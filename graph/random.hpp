// (C) 2014 Arek Olek

#pragma once

#include <random>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random_spanning_tree.hpp>
#include <boost/graph/named_function_params.hpp>

#include "range.hpp"

template <class Vertex, class Graph, class Generator>
Vertex random_neighbor(Vertex const & v, Graph const & G, Generator& gen) {
  uint d = out_degree(v, G)-1;
  auto k = std::uniform_int_distribution<uint>{0, d}(gen);
  auto it = adjacent_vertices(v, G).first;
  while(k--) ++it;
  return *it;
}

template <class Graph, class Tree>
Tree random_tree(Graph& G, unsigned seed) {
  std::default_random_engine gen(seed);
  unsigned n = num_vertices(G);
  std::vector<bool> visited(n, false);
  Tree T(n);
  unsigned v = std::uniform_int_distribution<unsigned>{0, n-1}(gen);
  while(num_edges(T) < n-1) {
    visited[v] = true;
    auto w = random_neighbor(v, G, gen);
    if(!visited[w]) add_edge(v, w, T);
    v = w;
  }
  return T;
}

template <class Graph, class Tree>
Tree wilson_tree(Graph const & G, unsigned seed) {
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex;
  unsigned n = num_vertices(G);
  std::vector<vertex> pred(n);
  std::default_random_engine gen(seed);
  random_spanning_tree(G, gen, boost::predecessor_map(&pred[0]));
  Tree T(n);
  for(unsigned i = 0; i < n; ++i)
    if(pred[i] != boost::graph_traits<Graph>::null_vertex())
      add_edge(i, pred[i], T);
  return T;
}
