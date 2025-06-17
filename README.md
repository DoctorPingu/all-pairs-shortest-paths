# C++ All-Pairs Shortest Paths Implementation

This project is a C++ implementation of algorithms designed to solve the All-Pairs Shortest Path (APSP) problem. Given a directed, weighted graph, these algorithms find the shortest path between every pair of vertices.

The repository includes two primary algorithms: Johnson's algorithm, which is highly efficient for sparse graphs, and the Floyd-Warshall algorithm, which is suitable for dense graphs. Both implementations are capable of handling negative edge weights and include robust logic for detecting negative-weight cycles, a condition under which shortest paths are not well-defined.

## Core Features

-   **Multiple Algorithms**: Provides both Johnson's algorithm and the Floyd-Warshall algorithm, allowing users to choose the most appropriate method based on graph density.
-   **Negative Weight Compatibility**: Correctly computes shortest paths in graphs that contain negative edge weights.
-   **Negative Cycle Detection**: Automatically detects the presence of negative-weight cycles. Johnson's algorithm will throw a `std::runtime_error` if a negative cycle is found.
-   **Generic Graph Representation**: The core `Graph` class is templated, enabling it to work with various numeric types such as `int` and `double`.
-   **Extensive Testing**: A comprehensive test suite using the Google Test framework is included to ensure the correctness and robustness of the implementations.
-   **File Parsing**: The `Graph` class can be constructed directly from a text file that defines the graph structure.

## Implementation Details

-   **Johnson's Algorithm**: This algorithm works by re-weighting the graph's edges to eliminate negative weights. This is achieved by computing a "potential" for each vertex using a Bellman-Ford-like procedure. With non-negative edges, Dijkstra's algorithm can then be efficiently run from each vertex to find the shortest paths.
-   **Floyd-Warshall Algorithm**: This algorithm uses a dynamic programming approach. It iteratively considers each vertex `k` and updates the shortest path between every pair of vertices `i` and `j` if the path from `i` to `j` via `k` is shorter.
-   **Graph Structure**: The graph is represented using an adjacency list (`std::vector<std::unordered_map<int, T>>`), which provides efficient access to the neighbors of any given vertex, making it ideal for sparse graphs.

## How to Build and Run Tests

The project comes with a full test suite built on the Google Test framework to validate the functionality.

### Prerequisites

-   A C++ compiler that supports the C++17 standard (e.g., GCC 9+, Clang 10+).
-   [Google Test](https://github.com/google/googletest) (gtest) installed and accessible by the compiler.

### Compiling and Running

1.  Navigate to the project's root directory in your terminal.
2.  Compile the source files using the following command. You may need to adjust the include (`-I`) and library (`-L`) paths to match your Google Test installation.

    ```bash
    g++ -std=c++17 main.cpp -o run_tests -I/path/to/googletest/include -L/path/to/googletest/lib -lgtest -lgtest_main -pthread
    ```

3.  Execute the compiled test runner:

    ```bash
    ./run_tests
    ```

    All tests should pass, confirming that the implementation is working correctly.

## Example Usage

Here is a simple example demonstrating how to use the `Graph` class and the APSP functions with the provided `tinyEWD.txt` dataset.

```cpp
#include "graph.hpp"
#include <iostream>
#include <vector>
#include <stdexcept>

int main() {
    // 1. Create a graph instance from an input file
    Graph<int> G("tinyEWD.txt");

    try {
        // 2. Calculate all-pairs shortest paths using Johnson's algorithm
        std::vector<std::vector<int>> distances = johnsonAPSP(G);

        // 3. Print a specific shortest path distance (e.g., from vertex 0 to 6)
        std::cout << "Shortest distance from vertex 0 to 6 is: " << distances << std::endl;

        // 4. Print another path (e.g., from vertex 3 to 5)
        std::cout << "Shortest distance from vertex 3 to 5 is: " << distances << std::endl;

    } catch (const std::runtime_error& e) {
        std::cerr << "Execution failed: " << e.what() << std::endl;
    }

    return 0;
}
```

### Expected Output
```cpp
Shortest distance from vertex 0 to 6 is: 151
Shortest distance from vertex 3 to 5 is: 154
```
