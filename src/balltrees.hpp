#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <limits>
#include <algorithm>
#include <iomanip>
#include<fstream>
#include<chrono>
#include "utils.hpp"

using namespace std;

// Class to represent a Ball (hypersphere) containing points
class Ball {
private:
    Point center;
    double radius;

public:
    Ball() : radius(0.0) {}
    
    Ball(const Point& c, double r) : center(c), radius(r) {}
    
    const Point& getCenter() const { return center; }
    
    double getRadius() const { return radius; }
    
    // Check if a point is inside the ball
    bool contains(const Point& p) const {
        return center.distanceTo(p) <= radius;
    }
    
    // Get minimum distance from a point to the ball surface
    double minDistanceTo(const Point& p) const {
        double dist = center.distanceTo(p);
        return max(0.0, dist - radius);
    }
};

// Node structure for Ball Tree
class BallTreeNode {
private:
    Ball ball;
    vector<Point> points;
    BallTreeNode* left;
    BallTreeNode* right;
    bool isLeaf;

public:
    BallTreeNode() : left(nullptr), right(nullptr), isLeaf(true) {}
    
    ~BallTreeNode() {
        delete left;
        delete right;
    }
    
    void setBall(const Ball& b) { ball = b; }
    
    const Ball& getBall() const { return ball; }
    
    void setPoints(const vector<Point>& pts) {
        points = pts;
        isLeaf = true;
    }
    
    const vector<Point>& getPoints() const { return points; }
    
    void setLeft(BallTreeNode* l) { left = l; isLeaf = false; }
    
    void setRight(BallTreeNode* r) { right = r; isLeaf = false; }
    
    BallTreeNode* getLeft() const { return left; }
    
    BallTreeNode* getRight() const { return right; }
    
    bool getIsLeaf() const { return isLeaf; }
};

// Main Ball Tree class
class BallTree {
private:
    BallTreeNode* root;
    int leafSize;
    int dimension;
    
    // Calculate centroid of a set of points
    Point calculateCentroid(const vector<Point>& pts) {
        if (pts.empty()) return Point(dimension);
        
        Point centroid(dimension);
        for (const auto& p : pts) {
            for (int i = 0; i < dimension; i++) {
                centroid[i] += p[i];
            }
        }
        
        for (int i = 0; i < dimension; i++) {
            centroid[i] /= pts.size();
        }
        
        return centroid;
    }
    
    // Create a ball containing all points
    Ball createBall(const vector<Point>& pts) {
        Point center = calculateCentroid(pts);
        double maxDist = 0.0;
        
        for (const auto& p : pts) {
            double dist = center.distanceTo(p);
            maxDist = max(maxDist, dist);
        }
        
        return Ball(center, maxDist);
    }
    
    // Find the dimension with maximum spread
    int findMaxSpreadDimension(const vector<Point>& pts) {
        int maxDim = 0;
        double maxSpread = 0.0;
        
        for (int d = 0; d < dimension; d++) {
            double minVal = numeric_limits<double>::max();
            double maxVal = numeric_limits<double>::lowest();
            
            for (const auto& p : pts) {
                minVal = min(minVal, p[d]);
                maxVal = max(maxVal, p[d]);
            }
            
            double spread = maxVal - minVal;
            if (spread > maxSpread) {
                maxSpread = spread;
                maxDim = d;
            }
        }
        
        return maxDim;
    }
    
    // Recursive function to build the ball tree
    BallTreeNode* buildTree(vector<Point>& pts) {
        if (pts.empty()) return nullptr;
        
        BallTreeNode* node = new BallTreeNode();
        Ball ball = createBall(pts);
        node->setBall(ball);
        
        // If points <= leafSize, make it a leaf node
        if (pts.size() <= leafSize) {
            node->setPoints(pts);
            return node;
        }
        
        // Find dimension with maximum spread
        int splitDim = findMaxSpreadDimension(pts);
        
        // Sort points along the split dimension
        sort(pts.begin(), pts.end(), [splitDim](const Point& a, const Point& b) {
            return a[splitDim] < b[splitDim];
        });
        
        // Split points into two groups
        int mid = pts.size() / 2;
        vector<Point> leftPoints(pts.begin(), pts.begin() + mid);
        vector<Point> rightPoints(pts.begin() + mid, pts.end());
        
        // Recursively build left and right subtrees
        node->setLeft(buildTree(leftPoints));
        node->setRight(buildTree(rightPoints));
        
        return node;
    }
    
    // Recursive function to print the tree
    void printTree(BallTreeNode* node, int depth) {
        if (!node) return;
        
        string indent(depth * 2, ' ');
        cout << indent << "Depth " << depth << ": ";
        cout << "Center = ";
        node->getBall().getCenter().print();
        cout << ", Radius = " << fixed << setprecision(3) 
             << node->getBall().getRadius();
        
        if (node->getIsLeaf()) {
            cout << " [LEAF with " << node->getPoints().size() << " points]" << endl;
        } else {
            cout << " [INTERNAL]" << endl;
            printTree(node->getLeft(), depth + 1);
            printTree(node->getRight(), depth + 1);
        }
    }
    
    // Recursive nearest neighbor search
    void searchNearest(BallTreeNode* node, const Point& query, 
                      Point& nearest, double& minDist) {
        if (!node) return;
        
        // Check if this ball can contain a closer point
        double ballDist = node->getBall().minDistanceTo(query);
        if (ballDist >= minDist) return;
        
        if (node->getIsLeaf()) {
            // Check all points in the leaf
            for (const auto& p : node->getPoints()) {
                double dist = query.distanceTo(p);
                if (dist < minDist) {
                    minDist = dist;
                    nearest = p;
                }
            }
        } else {
            // Search children, prioritizing the closer one
            BallTreeNode* first = node->getLeft();
            BallTreeNode* second = node->getRight();
            
            double leftDist = first ? first->getBall().minDistanceTo(query) 
                                    : numeric_limits<double>::max();
            double rightDist = second ? second->getBall().minDistanceTo(query) 
                                      : numeric_limits<double>::max();
            
            if (rightDist < leftDist) {
                swap(first, second);
            }
            
            searchNearest(first, query, nearest, minDist);
            searchNearest(second, query, nearest, minDist);
        }
    }
    
public:
    BallTree(int dim, int leaf = 1) : root(nullptr), dimension(dim), leafSize(leaf) {}
    
    ~BallTree() {
        delete root;
    }
    
    // Build the ball tree from a set of points
    void build(vector<Point>& points) {
        if (points.empty()) {
            cout << "Error: No points to build tree!" << endl;
            return;
        }
        
        root = buildTree(points);
        // cout << "\nBall Tree built successfully with " << points.size() 
            //  << " points!" << endl;
    }
    
    // Print the entire tree structure
    void print() {
        if (!root) {
            cout << "Tree is empty!" << endl;
            return;
        }
        
        cout << "\n========== Ball Tree Structure ==========" << endl;
        printTree(root, 0);
        cout << "=========================================" << endl;
    }
    
    // Find nearest neighbor to query point
    Point findNearestNeighbor(const Point& query) {
        if (!root) {
            cout << "Error: Tree is empty!" << endl;
            return Point(dimension);
        }
        
        Point nearest(dimension);
        double minDist = numeric_limits<double>::max();
        
        searchNearest(root, query, nearest, minDist);
        
        return nearest;
    }
    
    // Find distance to nearest neighbor
    double findDistanceToNearest(const Point& query) {
        Point nearest = findNearestNeighbor(query);
        return query.distanceTo(nearest);
    }
};

// Utility class to generate random points
// class PointGenerator {
// private:
//     mt19937 rng;
//     uniform_real_distribution<double> dist;

// public:
//     PointGenerator(double minVal = 0.0, double maxVal = 100.0) 
//         : rng(random_device{}()), dist(minVal, maxVal) {}
    
//     vector<Point> generate(int n, int k) {
//         vector<Point> points;
        
//         for (int i = 0; i < n; i++) {
//             Point p(k, i);
//             for (int j = 0; j < k; j++) {
//                 p[j] = dist(rng);
//             }
//             points.push_back(p);
//         }
        
//         return points;
//     }
    
//     Point generateSingle(int k) {
//         Point p(k);
//         for (int j = 0; j < k; j++) {
//             p[j] = dist(rng);
//         }
//         return p;
//     }
// };


// helper finctions
vector<Point> generateTestData(int n, int k, vector<int> &rarr ,int idx ,double minVal = -10.0, double maxVal = 10.0) {
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

void writeDataToCSV(int n, int k, int time ,const string& filename ) {
    ifstream inFile(filename);
    bool fileExists = inFile.good();
    inFile.close();
    
    ofstream outFile(filename, ios::app);  // Open in append mode
    
    if (!fileExists) {
        // Write header only if file is newly created
        outFile << "Number_of_Points,Dimensions,Time\n";
    }
    
    outFile << n << "," << k << "," << time << "\n";
    outFile.close();
}
