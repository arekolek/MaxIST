// (C) 2014 Arek Olek

#pragma once

#include <algorithm>
#include <stack>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "graph.hpp"
#include "range.hpp"

template <class Graph, class Tree>
void visit(unsigned v, Graph& G,
    std::vector<bool>& visited, std::vector<unsigned>& degree,
    Tree& T) {
  auto rdfs_rule = [&]() {
    unsigned best = 0, min = UINT_MAX;
    for(auto u : range(adjacent_vertices(v, G)))
      if(!visited[u] && degree[u] < min) {
        best = u;
        min = degree[u];
      }
    return best;
  };

  visited[v] = true;
  for(auto u : range(adjacent_vertices(v, G))) --degree[u];
  while(degree[v] > 0) {
    auto w = rdfs_rule();
    boost::add_edge(v, w, T);
    visit(w, G, visited, degree, T);
  }
}

template <class Graph, class Tree>
Tree rdfs_tree(Graph& G) {
  unsigned n = num_vertices(G);
  Tree T(n);
  std::vector<bool> visited(n, false);
  std::vector<unsigned> deg(n, 0);

  for(auto v : range(vertices(G))) deg[v] = out_degree(v, G);

  auto s = std::distance(deg.begin(), std::min_element(deg.begin(), deg.end()));
  visit(s, G, visited, deg, T);

  return T;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::default_random_engine g;
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
std::vector<unsigned> minima(Iter start, Iter end) {
  std::vector<unsigned> mins;
  auto min = *std::min_element(start, end);
  for(auto it = start; it != end; ++it)
    if(*it == min) mins.push_back(std::distance(start, it));
  return mins;
}

template <class Graph, class Tree>
void visit_node(unsigned v, Graph& G,
    std::vector<bool>& visited, std::vector<int>& degree,
    Tree& T) {
  auto rdfs_rule = [&]() {
    std::vector<unsigned> neighbors;
    int min = INT_MAX;
    for(auto u : range(adjacent_vertices(v, G)))
      if(!visited[u] && degree[u] <= min) {
        if(degree[u] < min) neighbors.clear();
        neighbors.push_back(u);
        min = degree[u];
      }
    return *select_randomly(neighbors.begin(), neighbors.end());
  };

  visited[v] = true;
  for(auto u : range(adjacent_vertices(v, G))) --degree[u];
  while(degree[v] > 0) {
    auto w = rdfs_rule();
    boost::add_edge(v, w, T);
    visit_node(w, G, visited, degree, T);
  }
}

template <class Tree, class Graph>
Tree rdfs_rand_tree(Graph& G) {
  unsigned n = num_vertices(G);
  Tree T(n);
  std::vector<bool> visited(n, false);
  std::vector<int> deg(n, 0);

  for(auto v : range(vertices(G))) deg[v] = out_degree(v, G);

  auto mins = minima(deg.begin(), deg.end());
  auto s = *select_randomly(mins.begin(), mins.end());

  visit_node(s, G, visited, deg, T);

  return T;
}

template <class Graph, class Tree>
Tree rdfs_best_tree(Graph& G) {
  auto upper = upper_bound(G);
  auto best_T = rdfs_rand_tree<Tree>(G);
  auto best_n = num_internal(best_T);
  for(int i = 0; i < 50; ++i) {
    if(best_n == upper) break;
    auto T = rdfs_rand_tree<Tree>(G);
    auto n = num_internal(T);
    if(n > best_n) {
      best_T = T;
      best_n = n;
    }
  }
  return best_T;
}
