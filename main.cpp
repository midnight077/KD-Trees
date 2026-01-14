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
#include "./src/kd_trees_ANN.hpp"
#include "./src/hnsw.hpp"
#include<omp.h>

using namespace std;

struct ComparisonResult {
    int index;
    double lsh_dist;
    double kd_dist;
    double kdANN_dist;
    double hnsw_dist;
    string filename;
};

int main(){
    int n , k , nqp;    // nqp -> number of query points
    cout<<"Input number of query points - ";
    cin>>n;
    cout<<endl<<"Input dimensions of each point - ";
    cin>>k;
    cout<<endl<<"Input number of query points - ";
    cin>>nqp;
    cout<<endl;

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
        saveTestDataToFile(querypoint_filename, nqp  , k , query_point);
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

    KDTreeANN treeANN(k);
    treeANN.build(points);

    HNSW hnsw(k, 16, 200);  // dimensions, M, efConstruction
    hnsw.build(points);
    hnsw.printStats();

    auto startTime = chrono::high_resolution_clock::now();
    
    vector<ComparisonResult> results(nqp);
    #pragma omp parallel for schedule(dynamic)
    for( int i = 0 ; i < nqp ; i++){
        auto result = lsh.findNearestNeighbor(query_point[i]);
        Point nearest_lsh = lsh.getPoint(result.first);
        // nearest_lsh.print();
        Point nearest_kd = tree.findNearestNeighbor(query_point[i]);
        // nearest_kd.print();
        // kd-ANN
        Point nearest_kd_ANN = treeANN.findNearestNeighbor(query_point[i]);

        // HNSW
        Point hnsw_nearest = hnsw.findNearestNeighbor(query_point[i], 50);
        
        
        
        double lsh_dist_from_query = query_point[i].distance(nearest_lsh);
        double kd_dist_from_query = query_point[i].distance(nearest_kd);
        double kdANN_dist_from_query = query_point[i].distance(nearest_kd_ANN);
        double hnsw_dist_from_query = query_point[i].distance(hnsw_nearest);
        // string difference_filename =  "./difference.csv";
        
        results[i] = {i , lsh_dist_from_query , kd_dist_from_query , kdANN_dist_from_query, hnsw_dist_from_query};
        
    }
    
    for(int i = 0; i < nqp; i++){
        // results[i].lshPoint.print();
        // results[i].kdPoint.print();
        string difference_filename =  "./onlydifference3.csv";
        writeComparisonToCSV(difference_filename, results[i].index, results[i].lsh_dist, results[i].kd_dist, results[i].kdANN_dist, results[i].hnsw_dist);
    }
    
    auto stopTime = chrono::high_resolution_clock::now();
    auto durationBuildTime = chrono::duration_cast<chrono::milliseconds>(stopTime - startTime);
    cout<<endl<<durationBuildTime.count()<<endl;
    std::cout << endl;
    return 0;
}