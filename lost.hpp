// (C) 2014 Arek Olek

#include <unordered_map>

#include <boost/graph/adjacency_matrix.hpp>

#include "debug.hpp"

namespace detail {
  boost::adjacency_matrix<boost::undirectedS> typedef graph;

  template <class Input, class Output>
  void copy_edges(const Input& in, Output& out) {
    auto es = edges(in);
    for(auto e = es.first; e != es.second; ++e)
      add_edge(source(*e, in), target(*e, in), out);
  }
}

template <class Graph>
class leaf_info {
public:
  leaf_info(Graph const & T_) : T(T_) {
    update();
  }
  bool is_path() const {
    return L.size() == 2;
  }
  std::vector<unsigned> const & leaves() const {
    return L;
  }
  unsigned branching(unsigned x) const {
    return B.at(x);
  }
  unsigned branching_neighbor(unsigned x) const {
    return N.at(x);
  }
  void update() {
    L.clear();
    B.clear();
    N.clear();
    auto vs = vertices(T);
    for(auto vit = vs.first; vit != vs.second; ++vit)
      if(degree(*vit, T) == 1) {
        unsigned l = *vit;
        L.push_back(l);
        unsigned p = l;
        unsigned x = *adjacent_vertices(l, T).first;
        while(degree(x, T) == 2) {
          auto it = adjacent_vertices(x, T).first;
          unsigned tmp = x;
          x = p == *it ? *(++it) : *it;
          p = tmp;
        }
        if(degree(x, T) > 2) {
          B[l] = x;
          N[l] = p;
        }
      }
  }
private:
  Graph const& T;
  std::vector<unsigned> L;
  std::unordered_map<unsigned, unsigned> B, N;
};

template <class Graph>
Graph prieto(Graph& G) {
  detail::graph M(num_vertices(G));
  detail::copy_edges(G, M);
  auto T = dfs_tree(G);
  leaf_info<Graph> info(T);
  int i = 0;
  do {
    show("tree" + std::to_string(i++) + ".dot", M, T);
  } while(!info.is_path() && rule2(M, T, info));
  return T;
}

template <class Graph, class Tree, class LeafInfo>
bool rule2(Graph& G, Tree& T, LeafInfo& info) {
  for(auto x : info.leaves())
    for(auto y : info.leaves())
      if(edge(x, y, G).second) {
        add_edge(x, y, T);
        remove_edge(info.branching(x), info.branching_neighbor(x), T);
        info.update();
        return true;
      }
  return false;
}
