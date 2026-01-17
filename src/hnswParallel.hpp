#pragma once

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_set>
#include <random>
#include <cmath>
#include <algorithm>
#include <limits>
#include <atomic>
#include <omp.h>
#include "utils.hpp"

using namespace std;

class HNSWParallel {
private:
    // Graph structure: layers[layer][nodeId] = vector of neighbor ids
    vector<vector<vector<int>>> layers;
    
    // All points stored in the structure
    vector<Point> points;
    
    // Parameters
    int M;           // Max number of connections per node
    int M0;          // Max connections for layer 0 (usually 2*M)
    int efConstruction;  // Size of dynamic candidate list during construction
    int maxLevel;    // Maximum level in the hierarchy
    double mL;       // Level multiplier (1/ln(M))
    int dimensions;  // Dimensionality of points
    
    // Entry point (node with highest level)
    atomic<int> entryPoint;
    atomic<int> entryPointLevel;
    
    // Synchronization primitives using OpenMP locks
    mutable omp_lock_t graphLock;              // For graph structure modifications
    mutable vector<omp_lock_t> nodeLocks;      // Per-node locks for fine-grained locking
    mutable omp_lock_t rngLock;                // For random number generation
    
    // Random number generator
    mt19937 rng;
    uniform_real_distribution<double> uniformDist;
    
    // Generate random level for a new node (thread-safe)
    int getRandomLevel() {
        omp_set_lock(&rngLock);
        double r = uniformDist(rng);
        int level = static_cast<int>(-log(r) * mL);
        omp_unset_lock(&rngLock);
        return min(level, maxLevel);
    }
    
    // Thread-local random level generation for parallel build
    int getRandomLevelThreadLocal(mt19937& localRng) {
        double r = uniformDist(localRng);
        int level = static_cast<int>(-log(r) * mL);
        return min(level, maxLevel);
    }
    
    // Priority queue element for search
    struct Candidate {
        double distance;
        int id;
        
        bool operator<(const Candidate& other) const {
            return distance < other.distance;  // Min-heap by default
        }
        
        bool operator>(const Candidate& other) const {
            return distance > other.distance;  // For max-heap
        }
    };
    
    // Search layer for nearest neighbors (thread-safe read)
    vector<Candidate> searchLayer(const Point& query, int entryPointId, 
                                   int ef, int layer) const {
        unordered_set<int> visited;
        
        // Min-heap for candidates (closest first)
        priority_queue<Candidate, vector<Candidate>, greater<Candidate>> candidates;
        
        // Max-heap for results (furthest first for easy removal)
        priority_queue<Candidate> results;
        
        double entryDist = query.distance(points[entryPointId]);
        candidates.push({entryDist, entryPointId});
        results.push({entryDist, entryPointId});
        visited.insert(entryPointId);
        
        while (!candidates.empty()) {
            Candidate current = candidates.top();
            candidates.pop();
            
            // If current candidate is further than the furthest result, stop
            if (current.distance > results.top().distance && 
                static_cast<int>(results.size()) >= ef) {
                break;
            }
            
            // Explore neighbors (read with lock)
            vector<int> neighbors;
            {
                omp_set_lock(&graphLock);
                if (layer < static_cast<int>(layers.size()) && 
                    current.id < static_cast<int>(layers[layer].size())) {
                    neighbors = layers[layer][current.id];
                }
                omp_unset_lock(&graphLock);
            }
            
            for (int neighborId : neighbors) {
                if (visited.find(neighborId) == visited.end()) {
                    visited.insert(neighborId);
                    
                    double neighborDist = query.distance(points[neighborId]);
                    
                    if (static_cast<int>(results.size()) < ef || 
                        neighborDist < results.top().distance) {
                        candidates.push({neighborDist, neighborId});
                        results.push({neighborDist, neighborId});
                        
                        if (static_cast<int>(results.size()) > ef) {
                            results.pop();
                        }
                    }
                }
            }
        }
        
        // Convert results to vector
        vector<Candidate> resultVec;
        while (!results.empty()) {
            resultVec.push_back(results.top());
            results.pop();
        }
        
        // Sort by distance (ascending)
        sort(resultVec.begin(), resultVec.end());
        
        return resultVec;
    }
    
    // Select neighbors using simple heuristic
    vector<int> selectNeighbors(const Point& query, vector<Candidate>& candidates, 
                                 int maxConnections) {
        vector<int> neighbors;
        
        // Sort candidates by distance
        sort(candidates.begin(), candidates.end());
        
        for (const auto& candidate : candidates) {
            if (static_cast<int>(neighbors.size()) >= maxConnections) {
                break;
            }
            neighbors.push_back(candidate.id);
        }
        
        return neighbors;
    }
    
    // Add bidirectional connection (thread-safe)
    void addConnection(int fromId, int toId, int layer) {
        int maxConn = (layer == 0) ? M0 : M;
        
        // Lock both nodes (in consistent order to avoid deadlock)
        int firstId = min(fromId, toId);
        int secondId = max(fromId, toId);
        
        omp_set_lock(&nodeLocks[firstId]);
        omp_set_lock(&nodeLocks[secondId]);
        
        // Add connection from -> to
        if (find(layers[layer][fromId].begin(), layers[layer][fromId].end(), toId) 
            == layers[layer][fromId].end()) {
            layers[layer][fromId].push_back(toId);
            
            // Prune if too many connections
            if (static_cast<int>(layers[layer][fromId].size()) > maxConn) {
                pruneConnectionsUnsafe(fromId, layer, maxConn);
            }
        }
        
        // Add connection to -> from
        if (find(layers[layer][toId].begin(), layers[layer][toId].end(), fromId) 
            == layers[layer][toId].end()) {
            layers[layer][toId].push_back(fromId);
            
            // Prune if too many connections
            if (static_cast<int>(layers[layer][toId].size()) > maxConn) {
                pruneConnectionsUnsafe(toId, layer, maxConn);
            }
        }
        
        omp_unset_lock(&nodeLocks[secondId]);
        omp_unset_lock(&nodeLocks[firstId]);
    }
    
    // Prune connections (unsafe version - caller must hold lock)
    void pruneConnectionsUnsafe(int nodeId, int layer, int maxConn) {
        vector<Candidate> neighbors;
        
        for (int neighborId : layers[layer][nodeId]) {
            double dist = points[nodeId].distance(points[neighborId]);
            neighbors.push_back({dist, neighborId});
        }
        
        sort(neighbors.begin(), neighbors.end());
        
        layers[layer][nodeId].clear();
        for (int i = 0; i < min(maxConn, static_cast<int>(neighbors.size())); i++) {
            layers[layer][nodeId].push_back(neighbors[i].id);
        }
    }

public:
    // Constructor
    HNSWParallel(int dims, int m = 16, int efConst = 200) 
        : M(m), M0(2 * m), efConstruction(efConst), maxLevel(16),
          dimensions(dims), entryPoint(-1), entryPointLevel(-1),
          rng(42), uniformDist(0.0, 1.0) {
        mL = 1.0 / log(static_cast<double>(M));
        omp_init_lock(&graphLock);
        omp_init_lock(&rngLock);
    }
    
    // Destructor to clean up locks
    ~HNSWParallel() {
        omp_destroy_lock(&graphLock);
        omp_destroy_lock(&rngLock);
        for (size_t i = 0; i < nodeLocks.size(); i++) {
            omp_destroy_lock(&nodeLocks[i]);
        }
    }
    
    // Build the HNSW graph from a set of points (parallel version)
    void build(const vector<Point>& inputPoints) {
        if (inputPoints.empty()) {
            return;
        }
        
        points.clear();
        layers.clear();
        entryPoint = -1;
        entryPointLevel = -1;
        
        int n = inputPoints.size();
        
        // Pre-allocate points vector
        points.reserve(n);
        
        // Pre-compute levels for all points
        vector<int> levels(n);
        for (int i = 0; i < n; i++) {
            levels[i] = getRandomLevel();
        }
        
        // Find maximum level
        int maxLevelFound = *max_element(levels.begin(), levels.end());
        
        // Pre-allocate layers structure
        layers.resize(maxLevelFound + 1);
        for (int l = 0; l <= maxLevelFound; l++) {
            layers[l].resize(n);
        }
        
        // Initialize node locks
        // First destroy any existing locks
        for (size_t i = 0; i < nodeLocks.size(); i++) {
            omp_destroy_lock(&nodeLocks[i]);
        }
        nodeLocks.resize(n);
        for (int i = 0; i < n; i++) {
            omp_init_lock(&nodeLocks[i]);
        }
        
        // Copy all points first
        points = inputPoints;
        
        // Find the entry point (node with highest level)
        int maxLevelIdx = 0;
        for (int i = 0; i < n; i++) {
            if (levels[i] > levels[maxLevelIdx]) {
                maxLevelIdx = i;
            }
        }
        entryPoint = maxLevelIdx;
        entryPointLevel = levels[maxLevelIdx];
        
        // Build incrementally but with parallel neighbor search
        // First point is already set as entry
        // Add remaining points with parallelization in batches
        
        int batchSize = max(1, n / (omp_get_max_threads() * 4));
        
        for (int batchStart = 0; batchStart < n; batchStart += batchSize) {
            int batchEnd = min(batchStart + batchSize, n);
            
            // Process batch in parallel
            #pragma omp parallel for schedule(dynamic)
            for (int i = batchStart; i < batchEnd; i++) {
                addPointParallel(i, levels[i]);
            }
        }
        
        cout << "HNSWParallel built with " << points.size() << " points and " 
             << layers.size() << " layers using " << omp_get_max_threads() << " threads" << endl;
    }
    
    // Add a single point to the graph (parallel-safe version)
    void addPointParallel(int newId, int level) {
        // Skip if this is the entry point (already handled)
        if (newId == entryPoint.load()) {
            return;
        }
        
        int currentEntryPoint = entryPoint.load();
        int currentEntryLevel = entryPointLevel.load();
        
        int currentPoint = currentEntryPoint;
        
        // Search from top layer down to level+1
        for (int l = currentEntryLevel; l > level; l--) {
            if (l < static_cast<int>(layers.size())) {
                auto results = searchLayer(points[newId], currentPoint, 1, l);
                if (!results.empty()) {
                    currentPoint = results[0].id;
                }
            }
        }
        
        // Insert into layers from level down to 0
        for (int l = min(level, currentEntryLevel); l >= 0; l--) {
            auto candidates = searchLayer(points[newId], currentPoint, efConstruction, l);
            
            int maxConn = (l == 0) ? M0 : M;
            auto neighbors = selectNeighbors(points[newId], candidates, maxConn);
            
            // Connect new point to neighbors
            for (int neighborId : neighbors) {
                addConnection(newId, neighborId, l);
            }
            
            if (!candidates.empty()) {
                currentPoint = candidates[0].id;
            }
        }
        
        // Update entry point if new point has higher level (atomic)
        int expected = currentEntryLevel;
        while (level > expected) {
            if (entryPointLevel.compare_exchange_weak(expected, level)) {
                entryPoint.store(newId);
                break;
            }
        }
    }
    
    // Original sequential addPoint for compatibility
    void addPoint(const Point& point) {
        omp_set_lock(&graphLock);
        
        int newId = points.size();
        points.push_back(point);
        
        // Get random level for this point
        int level = getRandomLevel();
        
        // Ensure we have enough layers
        while (static_cast<int>(layers.size()) <= level) {
            layers.push_back(vector<vector<int>>());
        }
        
        // Ensure each layer has space for this node
        for (int l = 0; l <= level; l++) {
            while (static_cast<int>(layers[l].size()) <= newId) {
                layers[l].push_back(vector<int>());
            }
        }
        
        while (static_cast<int>(layers[0].size()) <= newId) {
            layers[0].push_back(vector<int>());
        }
        
        // Add node lock
        nodeLocks.push_back(omp_lock_t());
        omp_init_lock(&nodeLocks.back());
        
        omp_unset_lock(&graphLock);  // Release graph lock for search operations
        
        // If this is the first point
        if (entryPoint == -1) {
            entryPoint = newId;
            entryPointLevel = level;
            return;
        }
        
        addPointParallel(newId, level);
    }
    
    // Find the nearest neighbor to a query point (thread-safe)
    Point findNearestNeighbor(const Point& query, int ef = 50) const {
        if (points.empty()) {
            throw runtime_error("HNSW is empty!");
        }
        
        int currentPoint = entryPoint.load();
        int currentLevel = entryPointLevel.load();
        
        // Traverse from top layer to layer 1
        for (int l = currentLevel; l >= 1; l--) {
            if (l < static_cast<int>(layers.size())) {
                auto results = searchLayer(query, currentPoint, 1, l);
                if (!results.empty()) {
                    currentPoint = results[0].id;
                }
            }
        }
        
        // Search layer 0 with ef candidates
        auto results = searchLayer(query, currentPoint, ef, 0);
        
        if (results.empty()) {
            return points[entryPoint.load()];
        }
        
        return points[results[0].id];
    }
    
    // Parallel batch query - find nearest neighbors for multiple queries
    vector<Point> findNearestNeighborBatch(const vector<Point>& queries, int ef = 50) const {
        int numQueries = queries.size();
        vector<Point> results(numQueries);
        
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < numQueries; i++) {
            results[i] = findNearestNeighbor(queries[i], ef);
        }
        
        return results;
    }
    
    // Find k nearest neighbors (thread-safe)
    vector<Point> findKNearestNeighbors(const Point& query, int k, int ef = 50) const {
        if (points.empty()) {
            throw runtime_error("HNSW is empty!");
        }
        
        ef = max(ef, k);  // ef should be at least k
        
        int currentPoint = entryPoint.load();
        int currentLevel = entryPointLevel.load();
        
        // Traverse from top layer to layer 1
        for (int l = currentLevel; l >= 1; l--) {
            if (l < static_cast<int>(layers.size())) {
                auto results = searchLayer(query, currentPoint, 1, l);
                if (!results.empty()) {
                    currentPoint = results[0].id;
                }
            }
        }
        
        // Search layer 0
        auto results = searchLayer(query, currentPoint, ef, 0);
        
        vector<Point> neighbors;
        for (int i = 0; i < min(k, static_cast<int>(results.size())); i++) {
            neighbors.push_back(points[results[i].id]);
        }
        
        return neighbors;
    }
    
    // Parallel batch k-NN query
    vector<vector<Point>> findKNearestNeighborsBatch(const vector<Point>& queries, 
                                                      int k, int ef = 50) const {
        int numQueries = queries.size();
        vector<vector<Point>> results(numQueries);
        
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < numQueries; i++) {
            results[i] = findKNearestNeighbors(queries[i], k, ef);
        }
        
        return results;
    }
    
    // Get point by index
    const Point& getPoint(int idx) const {
        return points[idx];
    }
    
    // Get number of points
    int size() const {
        return points.size();
    }
    
    // Get number of layers
    int getNumLayers() const {
        return layers.size();
    }
    
    // Set number of threads for parallel operations
    void setNumThreads(int numThreads) {
        omp_set_num_threads(numThreads);
    }
    
    // Get current number of threads
    int getNumThreads() const {
        return omp_get_max_threads();
    }
    
    // Print statistics
    void printStats() const {
        cout << "=== HNSWParallel Statistics ===" << endl;
        cout << "Number of points: " << points.size() << endl;
        cout << "Number of layers: " << layers.size() << endl;
        cout << "Entry point: " << entryPoint.load() << endl;
        cout << "Entry point level: " << entryPointLevel.load() << endl;
        cout << "M: " << M << ", M0: " << M0 << endl;
        cout << "Max threads: " << omp_get_max_threads() << endl;
        
        for (size_t l = 0; l < layers.size(); l++) {
            int nodesInLayer = 0;
            int totalConnections = 0;
            for (const auto& neighbors : layers[l]) {
                if (!neighbors.empty()) {
                    nodesInLayer++;
                    totalConnections += neighbors.size();
                }
            }
            cout << "Layer " << l << ": " << nodesInLayer << " nodes, " 
                 << totalConnections << " connections" << endl;
        }
    }
};