#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <unordered_map>
#include <limits>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <chrono>
#include<fstream>
#include<string>
#include "utils.hpp"

using namespace std;

// Class representing a single hash function (random projection)
class HashFunction {
private:
    Point randomVector;  // Random projection vector
    double offset;       // Random offset for binning
    double width;        // Bin width

public:
    HashFunction(int dimension, double binWidth , int rand) : width(binWidth) {
        // Generate random projection vector with Gaussian distribution
        random_device rd;
        mt19937 gen(rand);
        normal_distribution<double> dist(0.0, 1.0);
        
        vector<double> coords(dimension);
        for (int i = 0; i < dimension; i++) {
            coords[i] = dist(gen);
        }
        randomVector = Point(coords);
        
        // Random offset in [0, width)
        uniform_real_distribution<double> offsetDist(0.0, width);
        offset = offsetDist(gen);
    }
    
    // Hash a point to an integer bucket
    int hash(const Point& p) const {
        double projection = p.dot(randomVector);
        return static_cast<int>(floor((projection + offset) / width));
    }
};

// Class representing a hash table with multiple hash functions
class HashTable {
private:
    vector<HashFunction> hashFunctions;
    unordered_map<string, vector<int>> buckets;  // bucket key -> point indices
    int numHashFunctions;

public:
    HashTable(int numHashes, int dimension, double binWidth , int rand) 
        : numHashFunctions(numHashes) {
        for (int i = 0; i < numHashes; i++) {
            hashFunctions.push_back(HashFunction(dimension, binWidth, (i*rand)+10 ));
        }
    }
    
    // Generate bucket key for a point
    string getBucketKey(const Point& p) const {
        stringstream ss;
        for (int i = 0; i < numHashFunctions; i++) {
            ss << hashFunctions[i].hash(p);
            if (i < numHashFunctions - 1) ss << "_";
        }
        return ss.str();
    }
    
    // Insert a point into the hash table
    void insert(const Point& p, int pointIndex) {
        string key = getBucketKey(p);
        buckets[key].push_back(pointIndex);
    }
    
    // Get all point indices in the same bucket as query point
    vector<int> query(const Point& p) const {
        string key = getBucketKey(p);
        auto it = buckets.find(key);
        if (it != buckets.end()) {
            return it->second;
        }
        return vector<int>();
    }
    
    // Print bucket information
    void printBuckets(const vector<Point>& points) const {
        std::cout << "\n--- Hash Table Buckets ---" << endl;
        int bucketNum = 1;
        for (const auto& bucket : buckets) {
            std::cout << "Bucket " << bucketNum++ << " [Key: " << bucket.first << "]: ";
            std::cout << "Points {";
            for (size_t i = 0; i < bucket.second.size(); i++) {
                std::cout << bucket.second[i];
                if (i < bucket.second.size() - 1) 
                    std::cout << ", ";
            }
            std::cout << "}" << endl;
        }
    }
};

// Main LSH class
class LSH {
private:
    vector<Point> dataPoints;
    vector<HashTable> hashTables;
    int numTables;
    int numHashFunctions;
    int dimension;
    double binWidth;

public:
    LSH(int L, int k, int dim, double w = 4.0) 
        : numTables(L), numHashFunctions(k), dimension(dim), binWidth(w) {
        // Create L hash tables
        for (int i = 0; i < L; i++) {
            hashTables.push_back(HashTable(k, dim, w, i));
        }
    }
    
    // Add data points to LSH structure
    void addPoints(const vector<Point>& points) {
        dataPoints = points;
        
        // Insert each point into all hash tables
        for (size_t i = 0; i < points.size(); i++) {
            for (int j = 0; j < numTables; j++) {
                hashTables[j].insert(points[i], i);
            }
        }
    }
    
    // Print all buckets from all hash tables
    void printAllBuckets() const {
        for (int i = 0; i < numTables; i++) {
            std::cout << "\n=== Hash Table " << (i + 1) << " ===" << endl;
            hashTables[i].printBuckets(dataPoints);
        }
    }
    
    // Find nearest neighbor for a query point
    pair<int, double> findNearestNeighbor(const Point& query) const {
        if (dataPoints.empty()) {
            return {-1, -1.0};
        }
        
        // Collect candidate points from all hash tables
        unordered_map<int, bool> candidates;
        
        for (int i = 0; i < numTables; i++) {
            vector<int> bucketPoints = hashTables[i].query(query);
            for (int idx : bucketPoints) {
                candidates[idx] = true;
            }
        }
        
        // If no candidates found in buckets, do linear search
        if (candidates.empty()) {
            std::cout << "\nNo candidates found in LSH buckets. Performing linear search..." << endl;
            return linearSearchNearest(query);
        }
        
        // Find nearest among candidates
        int nearestIdx = -1;
        double minDistance = numeric_limits<double>::max();
        
        // cout << "\nCandidates from LSH buckets: {";
        bool first = true;
        for (const auto& cand : candidates) {
            // if (!first) cout << ", ";
            // cout << cand.first;
            first = false;
            
            double dist = dataPoints[cand.first].distance(query);
            if (dist < minDistance) {
                minDistance = dist;
                nearestIdx = cand.first;
            }
        }
        std::cout << "}" << endl;
        
        return {nearestIdx, minDistance};
    }
    
    // Linear search for nearest neighbor (fallback)
    pair<int, double> linearSearchNearest(const Point& query) const {
        int nearestIdx = -1;
        double minDistance = numeric_limits<double>::max();
        
        for (size_t i = 0; i < dataPoints.size(); i++) {
            double dist = dataPoints[i].distance(query);
            if (dist < minDistance) {
                minDistance = dist;
                nearestIdx = i;
            }
        }
        
        return {nearestIdx, minDistance};
    }
    
    // Get data point by index
    const Point& getPoint(int index) const {
        return dataPoints[index];
    }
    
    int getNumPoints() const { return dataPoints.size(); }
};

// int main() {
//     cout << "=== Locality Sensitive Hashing (LSH) for Nearest Neighbor Search ===" << endl;
//     cout << "================================================================\n" << endl;
    
//     int n=1000000;
//     vector<int> garr;
//     for(int i=0;i<= 200;i++){
//         garr.push_back(i);
//     }

//     for(int itr = 0 ; itr<10;itr++){
//         int k=2;
//         for(int i =0; i<10; i++){
//             // cout << "\nGenerating " << n << " random " << endl;

//             int idx = (itr*10) + i;
//             vector<Point> points = generateTestData(n,k,garr,idx );

//             // vector<Point> points = generateRandomPoints(n, k);
            
//             // cout << "\nGenerated Points:" << endl;
//             // for (int i = 0; i < n; i++) {
//             //     cout << "Point " << i << ": ";
//             //     points[i].print();
//             //     cout << endl;
//             // }
            
//             // Step 3: Create LSH structure and build hash tables
//             int L = 5;  // Number of hash tables
//             int numHashes = 4;  // Number of hash functions per table
//             double binWidth = 5.0;  // Bin width for hashing
            
//             // cout << "\n--- Building LSH Structure ---" << endl;
//             // cout << "Number of hash tables (L): " << L << endl;
//             // cout << "Hash functions per table: " << numHashes << endl;
//             // cout << "Bin width: " << binWidth << endl;
            
//             LSH lsh(L, numHashes, k, binWidth);
//             auto startBuildTime = chrono::high_resolution_clock::now();
//             lsh.addPoints(points);
//             auto stopBuildTime = chrono::high_resolution_clock::now();
//             auto durationBuildTime = chrono::duration_cast<chrono::milliseconds>(stopBuildTime - startBuildTime);
//             cout<<durationBuildTime.count();
//             cout << "\nLSH built successfully!" << endl;

//             string build_file_name = "lsh_build_Lkw"+ to_string(L) + to_string(numHashes) + to_string(binWidth) + ".csv";
//             writeDataToCSV(n,k,L,numHashes,durationBuildTime.count(), build_file_name);
            
            
//             // Print bucket information
//             // lsh.printAllBuckets();
            
//             // Step 4: Input query point
//             // cout << "\n\n=== Query Point ===" << endl;
//             Point queryPoint = generateQueryPoint(k, garr, idx);
//             // cout << "Query point: ";
//             // queryPoint.print();
//             // cout << endl;
            
//             // Step 5: Find nearest neighbor
//             // cout << "\n--- Finding Nearest Neighbor using LSH ---" << endl;
            
//             auto st = chrono::high_resolution_clock::now();
//             auto result = lsh.findNearestNeighbor(queryPoint);
//             auto stopt = chrono::high_resolution_clock::now();
//             auto dt = chrono::duration_cast<chrono::microseconds>(stopt - st);
//             cout<<dt.count();
            
//             string find_file_name = "lsh_find_Lkw"+ to_string(L) + to_string(numHashes) + to_string(binWidth) + ".csv";
//             writeDataToCSV(n,k,L,numHashes,dt.count(),find_file_name);

//             if (result.first == -1) {
//                 cout << "Error: Could not find nearest neighbor!" << endl;
//                 continue;
//             }
            
//             // cout << "\nNearest Neighbor: Point " << result.first << " -> ";
//             // lsh.getPoint(result.first).print();
//             // cout << endl;
            
//             // // Step 6: Print distance
//             // cout << "\n--- Distance Calculation ---" << endl;
//             // cout << "Distance between query point and nearest neighbor: " 
//             //     << fixed << setprecision(4) << result.second << endl;
            
//             // // Verification with linear search
//             // cout << "\n--- Verification (Linear Search) ---" << endl;
//             // auto trueResult = lsh.linearSearchNearest(queryPoint);
//             // cout << "True Nearest Neighbor: Point " << trueResult.first << " -> ";
//             // lsh.getPoint(trueResult.first).print();
//             cout << endl;
//             // cout << "True Minimum Distance: " << fixed << setprecision(4) 
//             //     << trueResult.second << endl;
            
//             // if (result.first == trueResult.first) {
//             //     cout << "LSH found the correct nearest neighbor!" << endl;
//             // } else {
//             //     cout << "LSH found an approximate neighbor (may not be exact)" << endl;
//             // }
//             k=k*2;
//         }
//     }
    
//     return 0;
// }