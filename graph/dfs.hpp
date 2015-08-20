// (C) 2014 Arek Olek

#pragma once

#include <vector>

#include <boost/graph/adjacency_list.hpp>


template <class Graph>
void visit(Graph const & G, unsigned v, std::vector<bool>& discovered, Graph& T) {
  discovered[v] = true;
  for(auto w : range(adjacent_vertices(v, G))) {
    if(!discovered[w]) {
      add_edge(v, w, T);
      visit(G, w, discovered, T);
    }
  }
}


template <class Graph>
Graph dfs_tree(Graph const & G) {
  std::vector<bool> visited(num_vertices(G));
  Graph tree(num_vertices(G));
  visit(G, 0, visited, tree);
  return tree;
}


