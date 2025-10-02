#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
#include <random>
#include <chrono>
#include <fstream>
#include <typeinfo>

using namespace std;

class Point {
public:
    vector<double> coords;
    
    Point() {}
    
    Point(const vector<double>& c) : coords(c) {}
    
    double distance(const Point& other) const {
        double dist = 0;
        for (size_t i = 0; i < coords.size(); i++) {
            double diff = coords[i] - other.coords[i];
            dist += diff * diff;
        }
        return sqrt(dist);
    }
    
    void print() const {
        cout << "(";
        for (size_t i = 0; i < coords.size(); i++) {
            cout << fixed << setprecision(2) << coords[i];
            if (i < coords.size() - 1) cout << ", ";
        }
        cout << ")";
    }
};

class KDNode {
public:
    Point point;
    KDNode* left;
    KDNode* right;
    int axis;
    
    KDNode(const Point& p, int a) : point(p), left(nullptr), right(nullptr), axis(a) {}
};

class KDTree {
private:
    KDNode* root;
    int k; // dimensions
    
    KDNode* buildTree(vector<Point>& points, int depth) {
        if (points.empty()) return nullptr;
        
        int axis = depth % k;
        
        // Calculate mean along the current axis
        double sum = 0;
        for (const Point& p : points) {
            sum += p.coords[axis];
        }
        double mean = sum / points.size();
        
        // Find the point closest to the mean
        int bestIdx = 0;
        double minDist = abs(points[0].coords[axis] - mean);
        for (int i = 1; i < points.size(); i++) {
            double dist = abs(points[i].coords[axis] - mean);
            if (dist < minDist) {
                minDist = dist;
                bestIdx = i;
            }
        }
        
        // Create node with the point closest to mean
        KDNode* node = new KDNode(points[bestIdx], axis);
        
        // Split points into left and right based on mean value
        vector<Point> leftPoints, rightPoints;
        for (int i = 0; i < points.size(); i++) {
            if (i == bestIdx) continue; // Skip the selected node point
            
            if (points[i].coords[axis] <= mean) {
                leftPoints.push_back(points[i]);
            } else {
                rightPoints.push_back(points[i]);
            }
        }
        
        // Recursively build subtrees
        node->left = buildTree(leftPoints, depth + 1);
        node->right = buildTree(rightPoints, depth + 1);
        
        return node;
    }
    
    void nearestNeighborSearch(KDNode* node, const Point& target, Point& best, double& bestDist) {
        if (!node) return;
        
        // Calculate distance to current node
        double dist = node->point.distance(target);
        
        // Update best if current is better
        if (dist < bestDist) {
            bestDist = dist;
            best = node->point;
        }
        
        // Determine which side to search first
        int axis = node->axis;
        bool goLeft = target.coords[axis] < node->point.coords[axis];
        
        KDNode* nearSide = goLeft ? node->left : node->right;
        KDNode* farSide = goLeft ? node->right : node->left;
        
        // Search the near side first
        nearestNeighborSearch(nearSide, target, best, bestDist);
        
        // Check if we need to search the far side
        double axisDistance = abs(target.coords[axis] - node->point.coords[axis]);
        if (axisDistance < bestDist) {
            nearestNeighborSearch(farSide, target, best, bestDist);
        }
    }

    int calculateHeight(KDNode* node) {
        if (!node) return 0;
        
        int leftHeight = calculateHeight(node->left);
        int rightHeight = calculateHeight(node->right);
        
        return 1 + max(leftHeight, rightHeight);
    }
    
    void deleteTree(KDNode* node) {
        if (node) {
            deleteTree(node->left);
            deleteTree(node->right);
            delete node;
        }
    }
    
public:
    KDTree(int dimensions) : root(nullptr), k(dimensions) {}
    
    ~KDTree() {
        deleteTree(root);
    }
    
    void build(vector<Point> points) {
        root = buildTree(points, 0);
    }

    int getHeight() {
        return calculateHeight(root);
    }
    
    Point findNearestNeighbor(const Point& target) {
        if (!root) {
            throw runtime_error("Tree is empty!");
        }
        
        Point best = root->point;
        double bestDist = root->point.distance(target);
        
        nearestNeighborSearch(root, target, best, bestDist);
        
        return best;
    }
    
    void printInOrder(KDNode* node) {
        if (node) {
            printInOrder(node->left);
            cout << "Point: ";
            node->point.print();
            cout << " (axis: " << node->axis << ")" << endl;
            printInOrder(node->right);
        }
    }
    
    void printTree() {
        cout << "K-D Tree structure:" << endl;
        printInOrder(root);
    }
};

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

// Brute force nearest neighbor search function
Point bruteForceNearestNeighbor(const vector<Point>& points, const Point& query) {
    if (points.empty()) {
        throw runtime_error("No points available for search!");
    }
    
    Point nearest = points[0];
    double minDistance = query.distance(points[0]);
    
    for (size_t i = 1; i < points.size(); i++) {
        double currentDistance = query.distance(points[i]);
        if (currentDistance < minDistance) {
            minDistance = currentDistance;
            nearest = points[i];
        }
    }
    
    return nearest;
}

void writeDataToCSV(int n, int k, int height, int time ,const string& filename = "mean_based.csv") {
    ifstream inFile(filename);
    bool fileExists = inFile.good();
    inFile.close();
    
    ofstream outFile(filename, ios::app);  // Open in append mode
    
    if (!fileExists) {
        // Write header only if file is newly created
        outFile << "Number_of_Points,Dimensions,Time,Tree_Height\n";
    }
    
    outFile << n << "," << k << "," << time <<"," << height << "\n";
    outFile.close();
}

int main() {

    vector<int> garr;
    for(int i=0;i<= 200;i++){
        garr.push_back(i);
    }

    int n=1000000;

    for(int itr = 0 ; itr<10;itr++){
        int k=2;
        for(int i =0; i<10; i++){
            
            int idx = (itr*10) + i;
            vector<Point> points = generateTestData(n,k,garr,idx );
                        
            // Build k-d tree
            KDTree tree(k);
            

            auto startBuildTime = chrono::high_resolution_clock::now();
            tree.build(points);
            auto stopBuildTime = chrono::high_resolution_clock::now();
            auto durationBuildTime = chrono::duration_cast<chrono::milliseconds>(stopBuildTime - startBuildTime);
            cout<<durationBuildTime.count();
            cout << "\nK-D Tree built successfully!" << endl;
            
            // cout<<endl<<endl;
            // tree.printTree();
            // cout<<endl<<endl;

            int h = tree.getHeight();
            cout << "height of kd tree : "<<h<<endl<<endl;
        
            writeDataToCSV(n,k,h,durationBuildTime.count());
            
            Point queryPoint= generateQueryPoint(k, garr, idx);
        
            // cout<<endl<<endl;
            // Point brute_force = bruteForceNearestNeighbor(points , queryPoint);
            // cout<<"Brute force solution ---------------- ";
            // brute_force.print();
            // cout<<endl<<endl;
            
            // Find nearest neighbor
            try {
                auto st = chrono::high_resolution_clock::now();
                Point nearest = tree.findNearestNeighbor(queryPoint);
                auto stopt = chrono::high_resolution_clock::now();
                auto dt = chrono::duration_cast<chrono::nanoseconds>(stopt - st);
    
                writeDataToCSV(n,k,h,dt.count(),"mean_based_find_point.csv");
                
                // cout << "\n=== RESULT ===" << endl;
                // cout << "Query point: ";
                // queryPoint.print();
                // cout << endl;
                
                // cout << "Nearest neighbor ---------------------------------";
                // nearest.print();
                // cout << endl<<endl;
                // cout << "Distance: " << fixed << setprecision(4) 
                //      << queryPoint.distance(nearest) << endl;
            }
            catch (const exception& e) {
                cout << "Error: " << e.what() << endl;
            }
            k = k*2;
        }    
    }
    return 0;
}
