// (C) 2014 Arek Olek

#pragma once

#include <deque>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "range.hpp"

template <class Graph, class Tree>
Tree greedy_tree(Graph const & G) {
  typedef std::pair<int, int> edge;

  unsigned n = num_vertices(G);
  Tree T(n);

  std::vector<bool> visited(n, false);
  std::deque<edge> edges;
  int v = random<int>(0, n-1), x, y;

  while(true) {
    visited[v] = true;

    int i = 0;
    for(auto w : shuffled(adjacent_vertices(v, G))) if(!visited[w]) {
      if(++i == 1) x = v, y = w;
      else         edges.emplace_back (v, w);
    }

    if(i == 0) {
      std::random_shuffle(edges.begin(), edges.end());
      do {
        if(edges.empty()) return T;
        std::tie(x, y) = edges.front();
        edges.pop_front();
      } while(visited[y]);
    }

    add_edge(x, y, T);
    v = y;
  }
}
