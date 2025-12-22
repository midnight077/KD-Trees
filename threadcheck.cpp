#include <iostream>
#include <thread>

int main() {
    // This tells you how many concurrent threads your CPU can truly support
    unsigned int n = std::thread::hardware_concurrency();
    std::cout << "Optimal number of threads: " << n << std::endl;
    return 0;
}