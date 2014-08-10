// (C) 2014 Arek Olek

#pragma once

#include <stack>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "range.hpp"

namespace detail {

  unsigned typedef node;
  std::pair<node, node> typedef edge;

  template <class Graph>
  void visit(Graph& G, Graph& T, node v,
    node& next_rank, std::vector<node>& rank,
    std::vector<node>& deg, std::stack<edge>& edges) {

    rank[v] = next_rank++;

    node w, x, y;
    unsigned min_deg = UINT_MAX;

    for(auto u : range(adjacent_vertices(v, G))) {
      --deg[u];

      if(rank[u] == 0 && deg[u] < min_deg) {
        w = u;
        min_deg = deg[u];
      }
    }

    if(min_deg == UINT_MAX) {
      while(!edges.empty() && rank[edges.top().second] != 0)
        edges.pop();
      if(edges.empty())
        return;
      std::tie(x, y) = edges.top();
    } else std::tie(x, y) = edge(v, w);

    add_edge(x, y, T);

    for(auto u : range(adjacent_vertices(x, G)))
      if(u != y && rank[u] == 0)
        edges.emplace(x, u);

    visit(G, T, y, next_rank, rank, deg, edges);
  }

}

template <class Graph>
Graph rdfs_tree(Graph& G) {
  unsigned n = num_vertices(G);
  Graph T(n);
  unsigned next_rank = 1;
  std::vector<detail::node> rank(n, 0), deg(n);
  std::stack<detail::edge> edges;

  for(auto v : range(vertices(G)))
    deg[v] = degree(v, G);

  detail::visit(G, T, 0, next_rank, rank, deg, edges);

  return T;
}
