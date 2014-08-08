// (C) 2014 Arek Olek

#include <functional>
#include <unordered_map>

#include <boost/functional/hash.hpp>
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

namespace std {
  template<typename S, typename T>
  struct hash<pair<S, T>> {
    inline size_t operator()(const pair<S, T> & v) const {
      size_t seed = 0;
      boost::hash_combine(seed, v.first);
      boost::hash_combine(seed, v.second);
      return seed;
    }
  };
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
  unsigned branching(unsigned l) const {
    return B.at(l);
  }
  unsigned parent(unsigned x, unsigned l) const {
    return P.at(uintpair(l, x));
  }
  void update() {
    L.clear();
    B.clear();
    P.clear();
    auto vs = vertices(T);
    for(auto vit = vs.first; vit != vs.second; ++vit)
      if(degree(*vit, T) == 1) {
        L.push_back(*vit);
        traverse(*vit, T);
      }
  }
  void traverse(unsigned l, Graph const & T) {
    traverse(l, l, next(l, l, T).second, T);
  }
  void traverse(unsigned l, unsigned a, unsigned b, Graph const & T) {
    P[uintpair(l, b)] = a;
    while(degree(b, T) == 2) {
      std::tie(a, b) = next(a, b, T);
      P[uintpair(l, b)] = a;
    }
    if(degree(b, T) > 2) {
      if(B.count(l) == 0)
        B[l] = b;
      auto vs = adjacent_vertices(b, T);
      for(auto v = vs.first; v != vs.second; ++v)
        if(*v != a)
          traverse(l, b, *v, T);
    }
  }
  std::pair<unsigned, unsigned> next(unsigned a, unsigned b, Graph const & T) {
    auto it = adjacent_vertices(b, T).first;
    return make_pair(b, a == *it ? *(++it) : *it);
  }
private:
  Graph const& T;
  std::vector<unsigned> L;
  std::pair<unsigned, unsigned> typedef uintpair;
  std::unordered_map<unsigned, unsigned> B;
  std::unordered_map<uintpair, unsigned> P;
};

template <class Graph>
Graph prieto(Graph& G) {
  //detail::graph M(num_vertices(G));
  //detail::copy_edges(G, M);
  auto T = dfs_tree(G);
  leaf_info<Graph> info(T);
  //int i = 0;
  do {
    //show("tree" + std::to_string(i++) + ".dot", M, T);
  } while(!info.is_path() && rule2(G, T, info));
  return T;
}

template <class Graph, class Tree, class LeafInfo>
bool rule2(Graph& G, Tree& T, LeafInfo& info) {
  for(auto l1 : info.leaves())
    for(auto l2 : info.leaves())
      if(edge(l1, l2, G).second) {
        add_edge(l1, l2, T);
        auto b = info.branching(l1);
        remove_edge(b, info.parent(b, l1), T);
        info.update();
        return true;
      }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule3(Graph& G, Tree& T, LeafInfo& info) {
  for(auto l : info.leaves()) {
    auto treeNeighbor = *adjacent_vertices(l, T).first;
    auto vs = adjacent_vertices(l, G);
    for(auto x = vs.first; x != vs.second; ++x)
      if(*x != treeNeighbor) {
        auto xl = info.parent(*x, l);
        if(degree(xl, T) > 2) {
          add_edge(l, *x, T);
          remove_edge(*x, xl, T);
          info.update();
          return true;
        }
      }
  }
  return false;
}

template <class Graph>
Graph lost_light(Graph& G) {
  auto T = dfs_tree(G);
  leaf_info<Graph> info(T);
  std::function<bool(Graph&,Graph&,leaf_info<Graph>&)> typedef rule;
  std::vector<rule> rules {
    rule2<Graph,Graph,leaf_info<Graph>>,
    rule3<Graph,Graph,leaf_info<Graph>>,
  };
  bool applied = true;
  while(applied && !info.is_path()) {
    applied = false;
    for(auto rule : rules) {
      if(rule(G, T, info)) {
        applied = true;
        break;
      }
    }
  }
  return T;
}
