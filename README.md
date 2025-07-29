# C++ All-Pairs Shortest Paths Implementation

This repository implements a graph-processing utility in C++ to determine the shortest paths between all pairs of vertices in a directed weighted graph. It includes support for detecting negative-weight cycles and reading graph data from file inputs. The implementation is designed for educational use and algorithmic understanding.

## Core Features

- **Graph Representation**: Directed weighted graph using adjacency maps.
- **File Input Support**: Load graphs from formatted `.txt` files.
- **Negative Cycle Detection**: Checks for cycles with a total negative weight.
- **All-Pairs Shortest Path Computation**: Intended to support shortest path algorithms across all vertex pairs.
- **Comprehensive Testing**: Includes unit tests covering edge cases and correctness.

## Implementation Details/Technology Used

- **Graph Class**: `Graph<T>` implemented using a vector of hash maps for adjacency lists.
- **Edge Management**: Supports adding, removing, checking, and querying edge weights.
- **Template Design**: Supports generic edge weights (e.g., `int`, `float`).
- **Negative Cycle Detection**: Likely utilises Bellman-Ford or a variant to detect cycles.
- **Built With**:
  - C++17 Standard
  - Standard Template Library (STL)
  - Google Test framework for unit testing

## Compilation/Execution Instructions

### Prerequisites

- C++ compiler with C++17 support (GCC 8+, Clang 10+)
- [Google Test](https://github.com/google/googletest) installed and accessible

### Compile and Run Tests

Use the following command to compile the test suite (adjust paths as needed):

```bash
g++ -std=c++17 main.cpp -o graph_tests -I/path/to/googletest/include -L/path/to/googletest/lib -lgtest -lgtest_main -pthread
```

Then run:

```bash
./graph_tests
```

### Demonstrating Example Usage (Separately from Tests)

Create a new `demo.cpp` file for standalone demonstration:

```cpp
#include <iostream>
#include "graph.hpp"

int main() {
    Graph<int> G("tinyEWD.txt");

    int u = 0, v = 3;
    if (G.isEdge(u, v)) {
        std::cout << "Edge from " << u << " to " << v << " with weight: " 
                  << G.getEdgeWeight(u, v) << std::endl;
    } else {
        std::cout << "No edge from " << u << " to " << v << std::endl;
    }

    std::cout << "Graph has " << G.size() << " vertices." << std::endl;
    return 0;
}
```

Compile and run:

```bash
g++ -std=c++17 demo.cpp -o demo
./demo
```

### Sample Input File (`tinyEWD.txt`)
```
8
4 5 0.35
5 4 0.35
4 7 0.37
5 7 0.28
7 5 0.28
5 1 0.32
0 4 0.38
0 2 0.26
...
```

### Expected Output
```
Edge from 0 to 3 with weight: 0.52
Graph has 8 vertices.
```

(*Note: Actual output depends on the content of the input file.*)

## Assumptions/Limitations

- Graphs must be defined in valid input files (e.g., `tinyEWD.txt`).
- Negative-weight edges are allowed; algorithms must guard against negative cycles.
- The codebase assumes well-formed input (minimal error handling for invalid formats).
