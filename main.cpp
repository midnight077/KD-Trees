#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <limits>
#include <algorithm>
#include <iomanip>
#include<fstream>
#include<chrono>
#include "./src/utils.hpp"
// #include "./src/balltrees.hpp"
#include "./src/lsh_euclidian.hpp"
#include "./src/kd_trees.hpp"
using namespace std;

int main(){
    int n , k , nqp;    // nqp -> number of query points
    cin>>n>>k>>nqp;

    vector<int> garr;
    for(int i=0;i<= 10;i++){
        garr.push_back(i);
    }
    //generate and save data set for both
    string database_filename = "./dataset/DB/"+ to_string(n)+"N"+to_string(k)+"K.txt";
    vector<Point> points;
    if(!fileExists(database_filename)){
        points = generateTestData(n, k, garr, 1, -100 , 100);
        saveTestDataToFile(database_filename, n , k , points );
    }
    else
        points = loadTestDataFromFile(database_filename);
    
    string querypoint_filename = "./dataset/Query/"+ to_string(nqp)+"N"+to_string(k)+"K.txt";

    vector<Point> query_point;
    if(!fileExists(querypoint_filename)){
        query_point = generateTestData(nqp , k, garr, 2, -100 , 100);
        saveTestDataToFile(querypoint_filename, n , k , query_point);
    }
    else    
        query_point = loadTestDataFromFile(querypoint_filename);    

    int L = 5;  // Number of hash tables
    int numHashes = 4;  // Number of hash functions per table
    double binWidth = 5.0;  // Bin width for hashing
            
    LSH lsh(L, numHashes, k, binWidth);
    lsh.addPoints(points);


    KDTree tree(k);
    tree.build(points);
        
    // auto startBuildTime = chrono::high_resolution_clock::now();
    // auto stopBuildTime = chrono::high_resolution_clock::now();
    // auto durationBuildTime = chrono::duration_cast<chrono::milliseconds>(stopBuildTime - startBuildTime);
    // std::cout<<durationBuildTime.count();
    // std::cout << "\nLSH built successfully!" << endl;

    // string build_file_name = "lsh_build_Lkw"+ to_string(L) + to_string(numHashes) + to_string(binWidth) + ".csv";
    // writeDataToCSV(n,k,L,numHashes,durationBuildTime.count(), build_file_name);
            
            
    // Print bucket information
    // lsh.printAllBuckets();
    // Step 4: Input query point
    for( int i = 0 ; i < nqp ; i++){
        auto result = lsh.findNearestNeighbor(query_point[i]);
        Point nearest_lsh = lsh.getPoint(result.first);
        nearest_lsh.print();
        Point nearest_kd = tree.findNearestNeighbor(query_point[i]);
        nearest_kd.print();
        int difference = nearest_lsh.distance(nearest_kd);

        writeComparisonToCSV("./difference.csv", query_point[i], nearest_lsh, nearest_kd, difference, i == 0);
    }
    
            // if (result.first == -1) {
            //     std::cout << "Error: Could not find nearest neighbor!" << endl;
            //     continue;
            // }
            
            // cout << "\nNearest Neighbor: Point " << result.first << " -> ";
            // lsh.getPoint(result.first).print();
            // cout << endl;
            
            // // Step 6: Print distance
            // cout << "\n--- Distance Calculation ---" << endl;
            // cout << "Distance between query point and nearest neighbor: " 
            //     << fixed << setprecision(4) << result.second << endl;
            
            // // Verification with linear search
            // cout << "\n--- Verification (Linear Search) ---" << endl;
            // auto trueResult = lsh.linearSearchNearest(queryPoint);
            // cout << "True Nearest Neighbor: Point " << trueResult.first << " -> ";
            // lsh.getPoint(trueResult.first).print();
            std::cout << endl;
    return 0;
}