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

template <typename T>
bool existsNegativeCycle(const Graph<T>& G) {              

  int nVertices = G.size();                                
  std::vector<T> dist(nVertices);

  bool updated;
  for (int iter = 0; iter < nVertices - 1; ++iter) {
    updated = false;
    for (int u = 0; u < nVertices; ++u) {
      for (const auto& [v, weight] : G.neighbours(u)) {

        if (dist[u] + weight < dist[v]) {
          dist[v] = dist[u] + weight;
          updated = true;
        }
      }
    }

    if (!updated) break;
  }

  for (int u = 0; u < nVertices; ++u) {
    for (const auto& [v, weight] : G.neighbours(u)) {
      if (dist[u] + weight < dist[v]) {
        return true;
      }
    }
  }

  return false;
}

template <typename T>
std::vector<std::vector<T>> johnsonAPSP(const Graph<T>& G) {
  int nVertices = G.size();
  T INF = infinity<T>();

  std::vector<std::vector<T>> distanceMatrix(nVertices, std::vector<T>(nVertices, INF));
  for (int i = 0; i < nVertices; ++i) {
    distanceMatrix[i][i] = 0;
  }

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

  for (int u = 0; u < nVertices; ++u) {
    for (const auto& [v, weight] : G.neighbours(u)) {
      if (potential[u] + weight < potential[v]) {
        throw std::runtime_error("Negative weight cycle detected, cannot compute Johnson's APSP");
      }
    }
  }

  for (int source = 0; source < nVertices; ++source) {

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

      for (const auto& [v, weight] : G.neighbours(u)) {
        T adjustedWeight = weight + potential[u] - potential[v];
        if (d + adjustedWeight < dist[v]) {
          dist[v] = d + adjustedWeight;
          pq.emplace(dist[v], v);
        }
      }
    }

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

template <typename T>
std::vector<std::vector<T>> floydWarshallAPSP(const Graph<T>& G) {
  int nVertices = G.size();
  T INF = infinity<T>();

  std::vector<std::vector<T>> dist(nVertices, std::vector<T>(nVertices, INF));
  for (int i = 0; i < nVertices; ++i) {
    dist[i][i] = 0;
  }

  for (int u = 0; u < nVertices; ++u) {
    for (const auto& [v, weight] : G.neighbours(u)) {
      if (weight < dist[u][v]) {
        dist[u][v] = weight;
      }
    }
  }

  for (int k = 0; k < nVertices; ++k) {
    for (int i = 0; i <nVertices; ++i) {
      if (dist[i][k] == INF) continue;
      for (int j = 0; j < nVertices; ++j) {
        if (dist[k][j] == INF) continue;
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
