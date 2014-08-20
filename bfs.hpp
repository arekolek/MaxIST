// (C) 2014 Arek Olek

#pragma once

#include <tuple>
#include <vector>
#include <queue>

template <class Graph>
Graph bfs_tree(Graph const & G) {
  typedef std::pair<int, int> edge;
  Graph T;
  int parent, v = 0;
  std::queue<edge> Q;
  std::vector<bool> V(num_vertices(G));
  Q.emplace(-1, v);
  V[v] = true;
  while(!Q.empty()) {
    std::tie(parent, v) = Q.front();
    Q.pop();
    for(auto w : range(adjacent_vertices(v, G))) {
      if(!V[w]) {
        V[w] = true;
        add_edge(v, w, T);
        Q.emplace(v, w);
      }
    }
  }
  return T;
}
