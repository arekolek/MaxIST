// (C) 2014 Arek Olek

#pragma once

#include <stack>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

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

    auto vs = adjacent_vertices(v, G);
    for(auto vit = vs.first; vit != vs.second; ++vit) {
      --deg[*vit];

      if(rank[*vit] == 0 && deg[*vit] < min_deg) {
        w = *vit;
        min_deg = deg[*vit];
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

    vs = adjacent_vertices(x, G);
    for(auto vit = vs.first; vit != vs.second; ++vit)
      if(*vit != y && rank[*vit] == 0)
        edges.emplace(x, *vit);

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

  auto vs = vertices(G);
  for(auto vit = vs.first; vit != vs.second; ++vit)
    deg[*vit] = degree(*vit, G);

  detail::visit(G, T, 0, next_rank, rank, deg, edges);

  return T;
}
