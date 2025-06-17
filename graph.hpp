/****************************************************************************************************
* README
* ---------
* This task was completed independently.
* As part of the assignment requirements, no GenAI were used to directly solve the assignment, write logic or produce code output.
*
* Any usage of copying and pasting was purely for:
* - Rearranging sections of existing code.
* - Reverting back to earlier versions of my own work through the submissions tab.
* - Improving clarity, formatting or correcting previous logic.
*
* All debugging and problem solving followed academic integrity principles outlined by UTS.
*
* References:
* - Bellman Ford Algorithm:
*   https://en.wikipedia.org/wiki/Bellman%E2%80%93Ford_algorithm
*
* - Floyd Warshall Algorithm:
*   https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm
*
* - Johnson Algorithm:
*   https://en.wikipedia.org/wiki/Johnson%27s_algorithm
****************************************************************************************************/

#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include <iostream>
#include <fstream>
#include <utility>
#include <functional>
#include <vector>
#include <string>
#include <queue>
#include <unordered_map>
#include <limits>

template <typename T>
class Graph {
 private:
  std::vector<std::unordered_map<int, T> > adjList {};
  int numVertices {};

 public:
  explicit Graph(int N);

  explicit Graph(const std::string& filename);

  void addEdge(int i, int j, T weight);

  void removeEdge(int i, int j);

  bool isEdge(int i, int j) const;

  T getEdgeWeight(int i, int j) const;

  int size() const;

  const std::unordered_map<int, T>& neighbours(int a) const {
    return adjList.at(a);
  }
};

template <typename T>
Graph<T>::Graph(int N) : adjList(N), numVertices {N} {}

template <typename T>
Graph<T>::Graph(const std::string& inputFile) {
  std::ifstream infile {inputFile};
  if (!infile) {
    std::cerr << inputFile << " could not be opened\n";
    return;
  }

  infile >> numVertices;
  adjList.resize(numVertices);
  int i {};
  int j {};
  double weight {};

  while (infile >> i >> j >> weight) {
    addEdge(i, j, static_cast<T>(weight));
  }
}

template <typename T>
int Graph<T>::size() const {
  return numVertices;
}

template <typename T>
void Graph<T>::addEdge(int i, int j, T weight) {
  if (i < 0 or i >= numVertices or j < 0 or j >= numVertices) {
    throw std::out_of_range("invalid vertex number");
  }
  adjList[i].insert({j, weight});
}

template <typename T>
void Graph<T>::removeEdge(int i, int j) {
  if (i >= 0 && i < numVertices && j >= 0 && j < numVertices) {
    adjList[i].erase(j);
  }
}

template <typename T>
bool Graph<T>::isEdge(int i, int j) const {
  if (i >= 0 && i < numVertices && j >= 0 && j < numVertices) {
    return adjList.at(i).contains(j);
  }
  return false;
}

template <typename T>
T Graph<T>::getEdgeWeight(int i, int j) const {
  return adjList.at(i).at(j);
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Graph<T>& G) {
  for (int i = 0; i < G.size(); ++i) {
    out << i << ':';
    for (const auto& [neighbour, weight] : G.neighbours(i)) {
      out << " (" << i << ", " << neighbour << ")[" << weight << ']';
    }
    out << '\n';
  }
  return out;
}

template <typename T>
T infinity() {
  if (std::numeric_limits<T>::has_infinity) {
    return std::numeric_limits<T>::infinity();
  } else {
    return std::numeric_limits<T>::max();
  }
}

// This function, existsNegativeCycle, uses the Bellman Ford relaxation method to detect if any negative weight cycle exists within the graph.
template <typename T>
bool existsNegativeCycle(const Graph<T>& G) {              

  // We initialise all vertec distances to 0.
  // This is so that every vertex is considered a starting point.
  int nVertices = G.size();                                
  std::vector<T> dist(nVertices);

  // Relax all edges (nVertices - 1) times.
  bool updated;
  for (int iter = 0; iter < nVertices - 1; ++iter) {
    updated = false;
    for (int u = 0; u < nVertices; ++u) {
      for (const auto& [v, weight] : G.neighbours(u)) {

        // Relaxation - Update dist[v] if a shorter path from u is found.
        if (dist[u] + weight < dist[v]) {
          dist[v] = dist[u] + weight;
          updated = true;
        }
      }
    }

    // If no distance was updated in this iteration, we can stop early.
    if (!updated) break;
  }

  // Perform one more pass to check for further improvements.
  // If any distance can be updated now, there is a negative cycle.
  for (int u = 0; u < nVertices; ++u) {
    for (const auto& [v, weight] : G.neighbours(u)) {
      if (dist[u] + weight < dist[v]) {
        return true;
      }
    }
  }

  // No further improvement means no negative weight cycle exists.
  return false;
}


// Handles all pairs paths even with negative edge weights. Has two phases:
// 1.) Reweighting phase
// 2.) Dijstra phase 
template <typename T>
std::vector<std::vector<T>> johnsonAPSP(const Graph<T>& G) {
  int nVertices = G.size();
  T INF = infinity<T>();

  // Preparee result matrix and initialise with infinity for no path, 0 on diagonals.
  std::vector<std::vector<T>> distanceMatrix(nVertices, std::vector<T>(nVertices, INF));
  for (int i = 0; i < nVertices; ++i) {
    distanceMatrix[i][i] = 0;
  }

  // 1.) Compute vertex potentials using Bellman Ford from a virtual source.
  std::vector<T> potential(nVertices, 0);
  for (int iter = 0; iter < nVertices - 1; ++iter) {
    bool changed = false;
    for (int u = 0; u < nVertices; ++u) {
      for (const auto& [v, weight] : G.neighbours(u)) {
        if (potential[u] + weight < potential[v]) {
          potential[v] = potential[u] + weight;
          changed = true;
        }
      }
    }
    if (!changed) break;
  }

  // Check for negative weight cycle in one more iteration.
  for (int u = 0; u < nVertices; ++u) {
    for (const auto& [v, weight] : G.neighbours(u)) {
      if (potential[u] + weight < potential[v]) {
        throw std::runtime_error("Negative weight cycle detected, cannot compute Johnson's APSP");
      }
    }
  }

  // 2.) Perform Dijkstra from each vertex on the re-weighted graph.
  for (int source = 0; source < nVertices; ++source) {

    // Min heap priority queue for (distance, vertex)
    using NodeDist = std::pair<T, int>;
    std::priority_queue<NodeDist, std::vector<NodeDist>, std::greater<NodeDist>> pq;
    std::vector<T> dist(nVertices, INF);
    dist[source] = 0;
    pq.emplace(0, source);

    while (!pq.empty()) {
      T d = pq.top().first;
      int u = pq.top().second;
      pq.pop();
      if (d != dist[u]) continue;

      // Relaxation: check all neighbours of u
      for (const auto& [v, weight] : G.neighbours(u)) {
        T adjustedWeight = weight + potential[u] - potential[v];
        if (d + adjustedWeight < dist[v]) {
          dist[v] = d + adjustedWeight;
          pq.emplace(dist[v], v);
        }
      }
    }

    // Convert distances back to original weights for the source.
    for (int v = 0; v < nVertices; ++v) {
      if (dist[v] == INF) {
        distanceMatrix[source][v] = INF;
      } else {
        distanceMatrix[source][v] = dist[v] + potential[v] - potential[source];
      }
    }
  }

  return distanceMatrix;
}


// Implementation of Floyd Warshall uses classic dynamic programmiing approach.
template <typename T>
std::vector<std::vector<T>> floydWarshallAPSP(const Graph<T>& G) {
  int nVertices = G.size();
  T INF = infinity<T>();

  // Initialise distance matrix with infinity.
  std::vector<std::vector<T>> dist(nVertices, std::vector<T>(nVertices, INF));
  for (int i = 0; i < nVertices; ++i) {
    dist[i][i] = 0;
  }

  // Set initial distances based on direct edges in the graph.
  for (int u = 0; u < nVertices; ++u) {
    for (const auto& [v, weight] : G.neighbours(u)) {
      if (weight < dist[u][v]) {
        dist[u][v] = weight;
      }
    }
  }

  // Floyd Warshall main loop - Try each vertex k as an intermediate
  for (int k = 0; k < nVertices; ++k) {
    for (int i = 0; i <nVertices; ++i) {
      // Skip if i -> k is unreachable.
      if (dist[i][k] == INF) continue;
      for (int j = 0; j < nVertices; ++j) {
        // Skip if k -> j is unreachable.
        if (dist[k][j] == INF) continue;
        
        // If a shorter path exists via k, update the distance.
        T newDistance = dist[i][k] + dist[k][j];
        if (newDistance < dist[i][j]) {
          dist[i][j] = newDistance;
        }
      }
    }
  }

  return dist;
}

#endif      // GRAPH_HPP_
