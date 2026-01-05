#pragma once

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_set>
#include <random>
#include <cmath>
#include <algorithm>
#include <limits>
#include "utils.hpp"

using namespace std;

class HNSW {
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
    int entryPoint;
    int entryPointLevel;
    
    // Random number generator
    mt19937 rng;
    uniform_real_distribution<double> uniformDist;
    
    // Generate random level for a new node
    int getRandomLevel() {
        double r = uniformDist(rng);
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
    
    // Search layer for nearest neighbors
    vector<Candidate> searchLayer(const Point& query, int entryPointId, 
                                   int ef, int layer) {
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
            
            // Explore neighbors
            if (layer < static_cast<int>(layers.size()) && 
                current.id < static_cast<int>(layers[layer].size())) {
                for (int neighborId : layers[layer][current.id]) {
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
    
    // Add bidirectional connection
    void addConnection(int fromId, int toId, int layer) {
        int maxConn = (layer == 0) ? M0 : M;
        
        // Add connection from -> to
        if (find(layers[layer][fromId].begin(), layers[layer][fromId].end(), toId) 
            == layers[layer][fromId].end()) {
            layers[layer][fromId].push_back(toId);
            
            // Prune if too many connections
            if (static_cast<int>(layers[layer][fromId].size()) > maxConn) {
                pruneConnections(fromId, layer, maxConn);
            }
        }
        
        // Add connection to -> from
        if (find(layers[layer][toId].begin(), layers[layer][toId].end(), fromId) 
            == layers[layer][toId].end()) {
            layers[layer][toId].push_back(fromId);
            
            // Prune if too many connections
            if (static_cast<int>(layers[layer][toId].size()) > maxConn) {
                pruneConnections(toId, layer, maxConn);
            }
        }
    }
    
    // Prune connections to keep only the closest ones
    void pruneConnections(int nodeId, int layer, int maxConn) {
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
    HNSW(int dims, int m = 16, int efConst = 200) 
        : M(m), M0(2 * m), efConstruction(efConst), maxLevel(16),
          dimensions(dims), entryPoint(-1), entryPointLevel(-1),
          rng(42), uniformDist(0.0, 1.0) {
        mL = 1.0 / log(static_cast<double>(M));
    }
    
    // Build the HNSW graph from a set of points
    void build(const vector<Point>& inputPoints) {
        if (inputPoints.empty()) {
            return;
        }
        
        points.clear();
        layers.clear();
        entryPoint = -1;
        entryPointLevel = -1;
        
        // Add points one by one
        for (const auto& point : inputPoints) {
            addPoint(point);
        }
        
        cout << "HNSW built with " << points.size() << " points and " 
             << layers.size() << " layers" << endl;
    }
    
    // Add a single point to the graph
    void addPoint(const Point& point) {
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
        
        // Also ensure layer 0 has space (all points exist in layer 0)
        while (static_cast<int>(layers[0].size()) <= newId) {
            layers[0].push_back(vector<int>());
        }
        
        // If this is the first point
        if (entryPoint == -1) {
            entryPoint = newId;
            entryPointLevel = level;
            return;
        }
        
        int currentPoint = entryPoint;
        
        // Search from top layer down to level+1
        for (int l = entryPointLevel; l > level; l--) {
            if (l < static_cast<int>(layers.size())) {
                auto results = searchLayer(point, currentPoint, 1, l);
                if (!results.empty()) {
                    currentPoint = results[0].id;
                }
            }
        }
        
        // Insert into layers from level down to 0
        for (int l = min(level, entryPointLevel); l >= 0; l--) {
            auto candidates = searchLayer(point, currentPoint, efConstruction, l);
            
            int maxConn = (l == 0) ? M0 : M;
            auto neighbors = selectNeighbors(point, candidates, maxConn);
            
            // Connect new point to neighbors
            for (int neighborId : neighbors) {
                addConnection(newId, neighborId, l);
            }
            
            if (!candidates.empty()) {
                currentPoint = candidates[0].id;
            }
        }
        
        // Update entry point if new point has higher level
        if (level > entryPointLevel) {
            entryPoint = newId;
            entryPointLevel = level;
        }
    }
    
    // Find the nearest neighbor to a query point
    Point findNearestNeighbor(const Point& query, int ef = 50) {
        if (points.empty()) {
            throw runtime_error("HNSW is empty!");
        }
        
        int currentPoint = entryPoint;
        
        // Traverse from top layer to layer 1
        for (int l = entryPointLevel; l >= 1; l--) {
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
            return points[entryPoint];
        }
        
        return points[results[0].id];
    }
    
    // Find k nearest neighbors
    vector<Point> findKNearestNeighbors(const Point& query, int k, int ef = 50) {
        if (points.empty()) {
            throw runtime_error("HNSW is empty!");
        }
        
        ef = max(ef, k);  // ef should be at least k
        
        int currentPoint = entryPoint;
        
        // Traverse from top layer to layer 1
        for (int l = entryPointLevel; l >= 1; l--) {
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
    
    // Print statistics
    void printStats() const {
        cout << "=== HNSW Statistics ===" << endl;
        cout << "Number of points: " << points.size() << endl;
        cout << "Number of layers: " << layers.size() << endl;
        cout << "Entry point: " << entryPoint << endl;
        cout << "Entry point level: " << entryPointLevel << endl;
        cout << "M: " << M << ", M0: " << M0 << endl;
        
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