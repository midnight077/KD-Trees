#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <limits>
#include <algorithm>
#include <iomanip>
#include<fstream>
#include<chrono>

using namespace std;

class Point {
private:
int id;

public:
    vector<double> coordinates;
    Point() : id(-1) {}
    
    Point(int dim, int pointId = -1) : coordinates(dim, 0.0), id(pointId) {}
    
    Point(const vector<double>& coords, int pointId = -1) 
        : coordinates(coords), id(pointId) {}
    
    int getDimension() const { return coordinates.size(); }
    
    int getId() const { return id; }
    
    void setId(int pointId) { id = pointId; }
    
    double& operator[](int index) { return coordinates[index]; }
    
    const double& operator[](int index) const { return coordinates[index]; }
    
    // Calculate Euclidean distance to another point
    double distance(const Point& other) const {
        if (coordinates.size() != other.coordinates.size()) {
            throw invalid_argument("Points must have same dimension");
        }
        
        double sum = 0.0;
        for (size_t i = 0; i < coordinates.size(); i++) {
            double diff = coordinates[i] - other.coordinates[i];
            sum += diff * diff;
        }
        return sqrt(sum);
    }
    
    // Dot product with another point/vector
    double dot(const Point& other) const {
        if (coordinates.size() != other.coordinates.size()) {
            throw invalid_argument("Points must have same dimension");
        }
        
        double result = 0.0;
        for (size_t i = 0; i < coordinates.size(); i++) {
            result += coordinates[i] * other.coordinates[i];
        }
        return result;
    }
    
    void print() const {
        cout << "(";
        for (size_t i = 0; i < coordinates.size(); i++) {
            cout << fixed << setprecision(2) << coordinates[i];
            if (i < coordinates.size() - 1) cout << ", ";
        }
        cout << ")";
    }
    
    const vector<double>& getCoordinates() const {
        return coordinates;
    }
};

vector<Point> generateTestData(int n, int k, vector<int> &rarr ,int idx ,double minVal = -100, double maxVal = 100.0) {
    // vector<int> rarr = {1,4,7,9,10,23,34,67,89,21,54,69,75};
    vector<Point> points;
    random_device rd;
    mt19937 gen(rarr[idx]); // Fixed seed for reproducible results
    uniform_real_distribution<double> dis(minVal, maxVal);
    
    for (int i = 0; i < n; i++) {
        vector<double> coords(k);
        for (int j = 0; j < k; j++) {
            coords[j] = dis(gen);
        }
        points.push_back(Point(coords));
    }  
    return points;
}

// Function to generate a random query point
Point generateQueryPoint(int k,vector<int> & rarr, int idx, double minVal = -10.0, double maxVal = 10.0) {
    // vector<int> rarr = {2,5,11,10,93,6,23,26,89,35,43,51,54,65,69};
    random_device rd;
    mt19937 gen(rarr[idx+1]); // Different seed for query point
    uniform_real_distribution<double> dis(minVal, maxVal);
    
    vector<double> coords(k);
    for (int j = 0; j < k; j++) {
        coords[j] = dis(gen);
    }
    
    return Point(coords);
}

// Function to generate and save test data to a file
void saveTestDataToFile(const string& filename, int n, int k, vector<Point> & points ) {
    
    
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << " for writing!" << endl;
        return;
    }
    
    // Write metadata: number of points and dimensions
    file << n << " " << k << endl;
    
    // Write each point
    for (const auto& point : points) {
        for (int i = 0; i < k; i++) {
            file << fixed << setprecision(10) << point[i];
            if (i < k - 1) file << " ";
        }
        file << endl;
    }
    
    file.close();
    cout << "Generated and Saved " << n << " points to " << filename << endl;
}

// Function to load points from a file into a vector<Point>
vector<Point> loadTestDataFromFile(const string& filename) {
    vector<Point> points;
    
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << " for reading!" << endl;
        return points;
    }
    
    int n, k;
    file >> n >> k;
    
    for (int i = 0; i < n; i++) {
        vector<double> coords(k);
        for (int j = 0; j < k; j++) {
            file >> coords[j];
        }
        points.push_back(Point(coords, i));  // Assign point ID as index
    }
    
    file.close();
    cout << "Loaded " << points.size() << " points from " << filename << endl;
    
    return points;
}

bool fileExists(const string& filename) {
    ifstream file(filename);
    return file.good();
}

void writeComparisonToCSV(const string& filename, int index, double lsh_dist, double kd_dist) {

    ifstream inFile(filename);
    bool fileExist = inFile.good();
    inFile.close();
    
    ofstream file;
    
    // Open in append mode
    file.open(filename, ios::app);
    
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << " for writing!" << endl;
        return;
    }
    
    // Write header if needed
    if (!fileExist) {
        file << "index,LSH_Distance,KD_Distance,Difference" << endl;
    }
    
    // Helper lambda to format point as string
    // auto pointToString = [](const Point& p) -> string {
    //     string result = "(";
    //     for (size_t i = 0; i < p.coordinates.size(); i++) {
    //         result += to_string(p[i]);
    //         if (i < p.coordinates.size() - 1) result += ";";
    //     }
    //     result += ")";
    //     return result;
    // };
    double difference = lsh_dist - kd_dist;
    // Write data row
    file << index << ","
         << fixed << setprecision(3) << lsh_dist << ","
         << fixed << setprecision(3) << kd_dist << ","
         << fixed << setprecision(3) << difference << endl;
    
    file.close();
}