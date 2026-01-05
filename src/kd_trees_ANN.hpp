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
#include "utils.hpp"
#include "kd_trees.hpp"

using namespace std;

class KDTreeANN {
private:
    KDNode* root;
    int k; // dimensions
    
    KDNode* buildTree(vector<Point>& points, int depth) {
        if (points.empty()) return nullptr;
        
        int axis = depth % k;
        
        // Sort points along the current axis using median
        sort(points.begin(), points.end(), [axis](const Point& a, const Point& b) {
            return a.coordinates[axis] < b.coordinates[axis];
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
        bool goLeft = target.coordinates[axis] < node->point.coordinates[axis];
        
        KDNode* nearSide = goLeft ? node->left : node->right;
        KDNode* farSide = goLeft ? node->right : node->left;
        
        // Search the near side first
        nearestNeighborSearch(nearSide, target, best, bestDist);
        
        // Check if we need to search the far side
        // double axisDistance = abs(target.coordinates[axis] - node->point.coordinates[axis]);
        // if (axisDistance < bestDist) {
        //     nearestNeighborSearch(farSide, target, best, bestDist);
        // }
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
    KDTreeANN(int dimensions) : root(nullptr), k(dimensions) {}
    
    ~KDTreeANN() {
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
            std::cout << "Point: ";
            node->point.print();
            std::cout << " (axis: " << node->axis << ")" << endl;
            printInOrder(node->right);
        }
    }
    
    void printTree() {
        cout<<"K-D Tree structure:" << endl;
        printInOrder(root);
    }
};