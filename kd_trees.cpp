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
        
        // Sort points along the current axis using median
        sort(points.begin(), points.end(), [axis](const Point& a, const Point& b) {
            return a.coords[axis] < b.coords[axis];
        });
        
        // Find median
        int median = points.size() / 2;
        
        // Create node
        KDNode* node = new KDNode(points[median], axis);
        
        // Split points into left and right subtrees
        vector<Point> leftPoints(points.begin(), points.begin() + median);
        vector<Point> rightPoints(points.begin() + median + 1, points.end());
        
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

vector<Point> generateTestData(int n, int k, double minVal = -50.0, double maxVal = 50.0) {
    vector<Point> points;
    random_device rd;
    mt19937 gen(rd()); // Fixed seed for reproducible results
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
Point generateQueryPoint(int k, double minVal = -10.0, double maxVal = 10.0) {
    random_device rd;
    mt19937 gen(rd()); // Different seed for query point
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

void writeDataToCSV(int n, int k, int height, int time ,const string& filename ) {
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
    int n, k=2;
    
    cout << "Enter number of points (n): ";
    cin >> n;

    // cout << "Enter number of dimensions (k): ";
    // cin >> k;
    for(int i =0; i<10; i++){

        vector<Point> points = generateTestData(n,k);
        
        // cout << "Enter " << n << " points (each point should have " << k << " coordinates):" << endl;
        // for (int i = 0; i < n; i++) {
        //     vector<double> coords(k);
        //     cout << "Point " << i + 1 << ": ";
        //     for (int j = 0; j < k; j++) {
        //         cin >> coords[j];
        //     }
        //     points.push_back(Point(coords));
        // }
        
        // Build k-d tree
        KDTree tree(k);
    
        auto startBuildTime = chrono::high_resolution_clock::now();
        tree.build(points);
        auto stopBuildTime = chrono::high_resolution_clock::now();
        auto durationBuildTime = chrono::duration_cast<chrono::microseconds>(stopBuildTime - startBuildTime);
        cout<<durationBuildTime.count();
        cout << "\nK-D Tree built successfully!" << endl;
        
        // Optional: Print tree structure
        // char choice;
        // cout << "Do you want to see the tree structure? (y/n): ";
        // cin >> choice;
        // if (choice == 'y' || choice == 'Y') {
        //     tree.printTree();
        // }
    
    
        int h = tree.getHeight();
        cout << "height of kd tree : "<<h<<endl<<endl;
    
        writeDataToCSV(n,k,h,durationBuildTime.count(), "median_based.csv");
        
        // Get query point
        // vector<double> queryCoords(k);
        // cout << "\nEnter the query point (" << k << " coordinates): ";
        // for (int j = 0; j < k; j++) {
        //     cin >> queryCoords[j];
        // }
        Point queryPoint= generateQueryPoint(k);
    
        Point brute_force = bruteForceNearestNeighbor(points , queryPoint);
        cout<<"Brute force solution - ";
        brute_force.print();
        
        // Find nearest neighbor
        try {
            auto st = chrono::high_resolution_clock::now();
            Point nearest = tree.findNearestNeighbor(queryPoint);
            auto stopt = chrono::high_resolution_clock::now();
            auto dt = chrono::duration_cast<chrono::microseconds>(stopt - st);

            writeDataToCSV(n,k,h,dt.count(),"median_based_find_point.csv");
            
            cout << "\n=== RESULT ===" << endl;
            cout << "Query point: ";
            // queryPoint.print();
            cout << endl;
            
            cout << "Nearest neighbor: ";
            // nearest.print();
            cout << endl;
            cout << "Distance: " << fixed << setprecision(4) 
                 << queryPoint.distance(nearest) << endl;
        }
        catch (const exception& e) {
            cout << "Error: " << e.what() << endl;
        }
        k = k*2;
    }    
    return 0;
}

/*
Example usage:

Input:
Enter number of dimensions (k): 2
Enter number of points (n): 5
Enter 5 points (each point should have 2 coordinates):
Point 1: 2 3
Point 2: 5 4
Point 3: 9 6
Point 4: 4 7
Point 5: 8 1

Enter the query point (2 coordinates): 6 5

Output:
Query point: (6.00, 5.00)
Nearest neighbor: (5.00, 4.00)
Distance: 1.4142

The algorithm now works by:
1. Building a k-d tree by recursively splitting points using MEAN along alternating axes
2. For each level, it calculates the mean value along the current axis
3. Selects the point closest to the mean as the splitting node
4. Divides remaining points: left subtree gets points <= mean, right subtree gets points > mean
5. For nearest neighbor search, it traverses the tree intelligently, pruning branches
   that cannot contain a closer point than the current best
6. Uses Euclidean distance for measuring proximity between points

Key differences from median-based approach:
- Mean-based splitting can create unbalanced trees if data is not uniformly distributed
- Points are split based on mean value rather than exact middle position
- The splitting node is chosen as the point closest to the calculated mean
- This approach may lead to different tree structures and potentially different performance characteristics
*/