// (C) 2014 Arek Olek

#pragma once

#include <stack>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

using namespace std;

unsigned typedef node;
pair<node, node> typedef edge;

template <class Graph>
void visit(Graph& G, Graph& T, node v,
  node& next_rank, vector<node>& rank,
  vector<node>& deg, stack<edge>& edges) {

  cerr << "visit " << v << endl;

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

  cerr << 1 << endl;

  if(min_deg == UINT_MAX) {
    while(!edges.empty() && rank[edges.top().second] != 0)
      edges.pop();
    if(edges.empty())
      return;
    tie(x, y) = edges.top();
  } else tie(x, y) = edge(v, w);

  cerr << 2 << endl;

  add_edge(x, y, T);

  vs = adjacent_vertices(x, G);
  for(auto vit = vs.first; vit != vs.second; ++vit)
    if(*vit != y && rank[*vit] == 0)
      edges.emplace(x, *vit);

  cerr << 3 << endl;

  visit(G, T, y, next_rank, rank, deg, edges);
}

template <class Graph>
Graph rdfs_tree(Graph& G) {
  unsigned n = num_vertices(G);
  Graph T(n);
  unsigned next_rank = 1;
  vector<node> rank(n, 0), deg(n);
  stack<edge> edges;

  auto vs = vertices(G);
  for(auto vit = vs.first; vit != vs.second; ++vit)
    deg[*vit] = degree(*vit, G);

  visit(G, T, 0, next_rank, rank, deg, edges);

  return T;
}
