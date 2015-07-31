// (C) 2014 Arek Olek

#pragma once

#include <random>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "range.hpp"

template <class IntType = int, class Generator = std::default_random_engine>
IntType random(IntType a = 0, IntType b = std::numeric_limits<IntType>::max()) {
  static Generator generator;
  return std::uniform_int_distribution<IntType>{a, b}(generator);
}

template <class Vertex, class Graph>
Vertex random_neighbor(Vertex const & v, Graph const & G) {
  unsigned target = random<unsigned>(0, degree(v, G)-1);
  for(auto w : range(adjacent_vertices(v, G)))
    if(target-- == 0) return w;
  throw "Ouch.";
}

template <class Graph>
Graph random_tree(Graph& G) {
  unsigned n = num_vertices(G);
  std::vector<bool> visited(n, false);
  Graph T(n);
  unsigned v = random<unsigned>(0, n-1);
  while(num_edges(T) < n-1) {
    visited[v] = true;
    auto w = random_neighbor(v, G);
    if(!visited[w]) add_edge(v, w, T);
    v = w;
  }
  return T;
}
