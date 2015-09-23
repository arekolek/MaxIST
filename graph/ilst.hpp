// (C) 2014 Arek Olek

#pragma once

#include <vector>

#include <boost/graph/adjacency_list.hpp>

template <class Graph>
class Ilst {
  std::vector<bool> discovered;
  unsigned conflicted, leaves;

  void visit(Graph const & G, Graph & tree, unsigned v) {
    discovered[v] = true;
    for(auto w : range(adjacent_vertices(v, G))) {
      if(!discovered[w]) {
        add_edge(v, w, tree);
        visit(G, tree, w);
      }
    }
    if(out_degree(v, tree) == 1) {
      ++leaves;
      if(edge(0, v, G).second) {
        conflicted = v;
      }
    }
  }
  std::pair<unsigned, unsigned> branch_edge(unsigned l, Graph const & tree) const {
    unsigned a = l, b = l, tmp;
    do {
      auto it = adjacent_vertices(b, tree).first;
      tmp = a;
      a = b;
      b = tmp == *it ? *(++it) : *it;
    } while(out_degree(b, tree) == 2);
    return std::make_pair(a, b);
  }
public:
  Ilst() : conflicted(0), leaves(0) {}

  Graph traverse(Graph const & G) {
    discovered.resize(num_vertices(G));
    Graph tree(num_vertices(G));
    visit(G, tree, 0);
    if(leaves > 2 && out_degree(0, tree) == 1 && conflicted != 0) {
      auto e = branch_edge(conflicted, tree);
      add_edge(0, conflicted, tree);
      remove_edge(e.first, e.second, tree);
    }
    return tree;
  }
};

template <class Graph>
Graph ilst(Graph const & G) {
  return Ilst<Graph>().traverse(G);
}
