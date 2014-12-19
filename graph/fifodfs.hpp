// (C) 2014 Arek Olek

#pragma once

#include <vector>
#include <deque>

#include "range.hpp"

template <class Graph>
Graph fifo_dfs_tree(Graph const & G) {
  typedef std::pair<int, int> edge;
  std::vector<bool> V(num_vertices(G));
  Graph T(num_vertices(G));
  std::deque<edge> Q;
  int parent, v = 0;
  Q.emplace_front(-1, v);
  while(!Q.empty()) {
    std::tie(parent, v) = Q.front();
    Q.pop_front();
    if(!V[v]) {
      V[v] = true;
      if(parent >= 0) add_edge(parent, v, T);
      int i = 0;
      for(auto w : range(adjacent_vertices(v, G))) if(!V[w]) {
        if(++i == 1) Q.emplace_front(v, w);
        else         Q.emplace_back (v, w);
      }
    }
  }
  return T;
}
