#include <iostream>
#include <vector>
#include <chrono>
#include "./src/utils.hpp"
#include "./src/hnsw.hpp"
#include "./src/hnswParallel.hpp"
#include "./src/kd_trees.hpp"

using namespace std;

int main() {
    int n, k, nqp;
    cout << "Input number of data points - ";
    cin >> n;
    cout << endl << "Input dimensions of each point - ";
    cin >> k;
    cout << endl << "Input number of query points - ";
    cin >> nqp;
    cout << endl;

    vector<int> garr;
    for (int i = 0; i <= 10; i++) {
        garr.push_back(i);
    }

    // Generate or load dataset
    string database_filename = "./dataset/DB/" + to_string(n) + "N" + to_string(k) + "K.txt";
    vector<Point> points;
    if (!fileExists(database_filename)) {
        points = generateTestData(n, k, garr, 1, -100, 100);
        saveTestDataToFile(database_filename, n, k, points);
    } else {
        points = loadTestDataFromFile(database_filename);
    }

    // Generate or load query points
    string querypoint_filename = "./dataset/Query/" + to_string(nqp) + "N" + to_string(k) + "K.txt";
    vector<Point> query_points;
    if (!fileExists(querypoint_filename)) {
        query_points = generateTestData(nqp, k, garr, 2, -100, 100);
        saveTestDataToFile(querypoint_filename, nqp, k, query_points);
    } else {
        query_points = loadTestDataFromFile(querypoint_filename);
    }

    // Build HNSW
    cout << "Building HNSW Parallel..." << endl;
    auto startBuildP = chrono::high_resolution_clock::now();
    
    HNSWParallel hnswP(k, 16, 200);  // dimensions, M, efConstruction
    hnswP.build(points);
    
    auto stopBuildP = chrono::high_resolution_clock::now();
    auto buildTimeP = chrono::duration_cast<chrono::milliseconds>(stopBuildP - startBuildP);
    cout << "HNSW Parallel build time: " << buildTimeP.count() << " ms" << endl;

    cout << "Building HNSW..." << endl;
    auto startBuild = chrono::high_resolution_clock::now();
    
    HNSW hnsw(k, 16, 200);  // dimensions, M, efConstruction
    hnsw.build(points);
    
    auto stopBuild = chrono::high_resolution_clock::now();
    auto buildTime = chrono::duration_cast<chrono::milliseconds>(stopBuild - startBuild);
    cout << "HNSW build time: " << buildTime.count() << " ms" << endl;
    
    hnsw.printStats();

    // Build KD-Tree for comparison
    cout << "\nBuilding KD-Tree..." << endl;
    KDTree kdtree(k);
    kdtree.build(points);

    // Query and compare
    cout << "\nQuerying..." << endl;
    auto startQuery = chrono::high_resolution_clock::now();
    
    int correct = 0;
    double totalHNSWDist = 0, totalKDDist = 0;
    
    for (int i = 0; i < nqp; i++) {
        Point hnsw_nearest = hnsw.findNearestNeighbor(query_points[i], 50);
        Point kd_nearest = kdtree.findNearestNeighbor(query_points[i]);
        
        double hnsw_dist = query_points[i].distance(hnsw_nearest);
        double kd_dist = query_points[i].distance(kd_nearest);
        
        cout<<endl<<i<<" "<<hnsw_dist-kd_dist<<endl;    
        totalHNSWDist += hnsw_dist;
        totalKDDist += kd_dist;
        
        // Check if HNSW found the exact nearest neighbor
        if (abs(hnsw_dist - kd_dist) < 1e-9) {
            correct++;
        }
    }
    
    auto stopQuery = chrono::high_resolution_clock::now();
    auto queryTime = chrono::duration_cast<chrono::milliseconds>(stopQuery - startQuery);
    
    cout << "\n=== Results ===" << endl;
    cout << "Query time: " << queryTime.count() << " ms" << endl;
    cout << "Accuracy (exact matches): " << (100.0 * correct / nqp) << "%" << endl;
    cout << "Average HNSW distance: " << (totalHNSWDist / nqp) << endl;
    cout << "Average KD-Tree distance: " << (totalKDDist / nqp) << endl;
    cout << "Average distance ratio (HNSW/KD): " << (totalHNSWDist / totalKDDist) << endl;

    return 0;
}