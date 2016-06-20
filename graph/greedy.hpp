// (C) 2014 Arek Olek

#pragma once

#include <deque>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "range.hpp"

template <class Graph, class Tree>
Tree greedy_tree(Graph const & G, unsigned seed) {
  typedef std::pair<int, int> edge;

  unsigned n = num_vertices(G);
  Tree T(n);

  std::vector<bool> visited(n, false);
  std::deque<edge> edges;
  std::default_random_engine gen(seed);
  unsigned v = std::uniform_int_distribution<unsigned>{0, n-1}(gen);
  int x = -1, y = -1;

  while(true) {
    visited[v] = true;

    int i = 0;
    for(auto w : shuffled(adjacent_vertices(v, G), gen)) if(!visited[w]) {
      if(++i == 1) x = v, y = w;
      else         edges.emplace_back (v, w);
    }

    if(i == 0) {
      std::shuffle(edges.begin(), edges.end(), gen);
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
