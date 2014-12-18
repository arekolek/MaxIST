// (C) 2014 Arek Olek

#pragma once

#include <algorithm>
#include <stack>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "range.hpp"

namespace detail {

  unsigned typedef node;
  std::pair<node, node> typedef edge;

  template <class Graph>
  void visit(Graph& G, Graph& T, node v,
    node& next_rank, std::vector<unsigned>& rank,
    std::vector<unsigned>& deg, std::stack<edge>& edges) {

    rank[v] = next_rank++;

    node w = -1, x, y;
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
  std::vector<unsigned> rank(n, 0), deg(n, 0);
  std::stack<detail::edge> edges;

  for(auto v : range(vertices(G)))
    deg[v] = degree(v, G);

  detail::visit(G, T, 0, next_rank, rank, deg, edges);

  return T;
}

template <class Graph>
Graph test_tree(Graph const & G) {
  typedef std::pair<int, int> edge;
  std::vector<unsigned> deg(num_vertices(G));
  std::vector<bool> V(num_vertices(G));
  std::stack<edge> Q;
  Graph T(num_vertices(G));
  int parent, v = 0;
  Q.emplace(-1, v);
  for(auto v : range(vertices(G))) deg[v] = degree(v, G);
  while(!Q.empty()) {
    std::tie(parent, v) = Q.top();
    Q.pop();
    if(!V[v]) {
      V[v] = true;
      if(parent >= 0) add_edge(parent, v, T);
      std::vector<unsigned> neighbors;
      for(auto w : range(adjacent_vertices(v, G))) if(!V[w]) {
        --deg[w];
        neighbors.push_back(w);
      }
      std::sort(neighbors.begin(), neighbors.end(), [&deg](unsigned a, unsigned b) {
        return deg[a] > deg[b];
      });
      for(auto w : neighbors) Q.emplace(v, w);
    }
  }
  return T;
}
