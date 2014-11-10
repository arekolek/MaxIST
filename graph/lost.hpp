// (C) 2014 Arek Olek

#include <functional>
#include <unordered_map>

#include <boost/functional/hash.hpp>
#include <boost/graph/adjacency_matrix.hpp>

#include "debug.hpp"
#include "range.hpp"

namespace detail {
  boost::adjacency_matrix<boost::undirectedS> typedef graph;
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
  unsigned branching_neighbor(unsigned l, unsigned x) const {
    return BN.at(uintpair(l, x));
  }
  bool on_branch(unsigned l, unsigned x) const {
    return BN.count(uintpair(l, x)) == 0;
  }
  void update() {
    L.clear();
    B.clear();
    P.clear();
    BN.clear();
    for(auto v : range(vertices(T)))
      if(degree(v, T) == 1) {
        L.push_back(v);
        traverse(v, T);
      }
  }
  void traverse(unsigned l, Graph const & T) {
    traverse(l, l, l, T);
  }
  void traverse(unsigned l, unsigned a, unsigned b, Graph const & T) {
    do {
      std::tie(a, b) = next(a, b, T);
      P[uintpair(l, b)] = a;
    } while(degree(b, T) == 2);
    if(degree(b, T) > 2) {
      B[l] = b;
      for(auto v : range(adjacent_vertices(b, T)))
        if(v != a)
          traverse(l, v, b, v, T);
    }
  }
  void traverse(unsigned l, unsigned blx, unsigned a, unsigned b, Graph const & T) {
    P[uintpair(l, b)] = a;
    BN[uintpair(l, b)] = blx;
    while(degree(b, T) == 2) {
      std::tie(a, b) = next(a, b, T);
      P[uintpair(l, b)] = a;
      BN[uintpair(l, b)] = blx;
    }
    if(degree(b, T) > 2) {
      for(auto v : range(adjacent_vertices(b, T)))
        if(v != a)
          traverse(l, blx, b, v, T);
    }
  }
  std::pair<unsigned, unsigned> next(unsigned a, unsigned b, Graph const & T) {
    auto it = adjacent_vertices(b, T).first;
    return std::make_pair(b, a == *it ? *(++it) : *it);
  }
private:
  Graph const& T;
  std::vector<unsigned> L;
  std::pair<unsigned, unsigned> typedef uintpair;
  std::unordered_map<unsigned, unsigned> B;
  std::unordered_map<uintpair, unsigned> P, BN;
};

class prieto {
public:
  template<class Graph, class Tree>
  int operator()(Graph& G, Tree& T) {
    //detail::graph M(num_vertices(G));
    //detail::copy_edges(G, M);
    leaf_info<Tree> info(T);
    int i = -1;
    do {
      ++i;
      //show("tree" + std::to_string(i++) + ".dot", M, T);
    } while(!info.is_path() && rule2(G, T, info));
    return i;
  }
};

template <class Graph, class Tree, class LeafInfo>
bool rule1(Graph& G, Tree& T, LeafInfo& info) {
  for(auto l1 : info.leaves())
    for(auto l2 : info.leaves())
      if(edge(l1, l2, G).second) {
        auto x = l1;
        auto y = *adjacent_vertices(l1, T).first;
        while(y != l2) {
          x = y;
          y = info.parent(y, l2);
          if(degree(x, T) > 2 && degree(y, T) > 2) {
            add_edge(l1, l2, T);
            remove_edge(x, y, T);
            info.update();
            return true;
          }
        }
      }
  return false;
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
    for(auto x : range(adjacent_vertices(l, G)))
      if(x != treeNeighbor) {
        auto xl = info.parent(x, l);
        if(degree(xl, T) > 2) {
          add_edge(l, x, T);
          remove_edge(x, xl, T);
          info.update();
          return true;
        }
      }
  }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule4(Graph& G, Tree& T, LeafInfo& info) {
  int n = num_vertices(G);
  std::vector<std::vector<std::pair<int,int>>> extra(n);
  for(auto l : info.leaves()) {
    auto treeNeighbor = *adjacent_vertices(l, T).first;
    for(auto x : range(adjacent_vertices(l, G)))
      if(x != treeNeighbor) {
        auto xl = info.parent(x, l);
        if(degree(xl, T) == 2)
          extra[xl].emplace_back(l, x);
      }
  }

  for(int l2 : info.leaves())
    for(int xl = 0; xl < n; ++xl)
      if(!extra[xl].empty() && edge(l2, xl, G).second) {
        int l, x;
        if(extra[xl].size() == 1 && extra[xl].begin()->first == l2) continue;
        std::tie(l, x) = extra[xl].begin()->first == l2 ? *(extra[xl].begin()+1) : *extra[xl].begin();
        add_edge(l, x, T);
        remove_edge(x, xl, T);
        info.update();

        add_edge(l2, xl, T);
        auto b = info.branching(l2);
        remove_edge(b, info.parent(b, l2), T);
        info.update();
        return true;
      }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule5(Graph& G, Tree& T, LeafInfo& info) {
  for(auto l : info.leaves()) {
    auto treeNeighbor = *adjacent_vertices(l, T).first;
    for(auto x : range(adjacent_vertices(l, G)))
      if(x != treeNeighbor && !info.on_branch(l, x)) {
        auto bl = info.branching(l);
        auto blx = info.branching_neighbor(l, x);
        if(degree(blx, T) > 2) {
          add_edge(l, x, T);
          remove_edge(bl, blx, T);
          info.update();
          return true;
        }
      }
  }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule6(Graph& G, Tree& T, LeafInfo& info) {
  int n = num_vertices(G);
  std::vector<std::vector<std::pair<int,int>>> extra(n);
  for(auto l : info.leaves()) {
    auto treeNeighbor = *adjacent_vertices(l, T).first;
    for(auto x : range(adjacent_vertices(l, G)))
      if(x != treeNeighbor && !info.on_branch(l, x)) {
        auto blx = info.branching_neighbor(l, x);
        if(degree(blx, T) == 2) {
          extra[blx].emplace_back(l, x);
        }
      }
  }

  for(int l2 : info.leaves())
    for(int blx = 0; blx < n; ++blx)
      if(!extra[blx].empty() && edge(l2, blx, G).second) {
        int l, x;
        if(extra[blx].size() == 1 && extra[blx].begin()->first == l2) continue;
        std::tie(l, x) = extra[blx].begin()->first == l2 ? *(extra[blx].begin()+1) : *extra[blx].begin();
        add_edge(l, x, T);
        remove_edge(info.branching(l), blx, T);
        info.update();

        add_edge(l2, blx, T);
        auto b = info.branching(l2);
        remove_edge(b, info.parent(b, l2), T);
        info.update();
        return true;
      }
  return false;
}

class lost_light {
public:
  template <class Graph, class Tree>
  int operator() (Graph& G, Tree& T) {
    leaf_info<Tree> info(T);
    std::function<bool(Graph&,Tree&,leaf_info<Tree>&)> typedef rule;
    std::vector<rule> rules {
      rule2<Graph,Tree,leaf_info<Tree>>,
      rule3<Graph,Tree,leaf_info<Tree>>,
      rule4<Graph,Tree,leaf_info<Tree>>,
      rule5<Graph,Tree,leaf_info<Tree>>,
      rule6<Graph,Tree,leaf_info<Tree>>,
    };
    int i = 0;
    bool applied = true;
    while(applied && !info.is_path()) {
      applied = false;
      for(auto rule : rules) {
        if(rule(G, T, info)) {
          ++i;
          applied = true;
          break;
        }
      }
    }
    return i;
  }
};

class lost {
public:
  template <class Graph, class Tree>
  int operator() (Graph& G, Tree& T) {
    leaf_info<Tree> info(T);
    std::function<bool(Graph&,Tree&,leaf_info<Tree>&)> typedef rule;
    std::vector<rule> rules {
      rule1<Graph,Tree,leaf_info<Tree>>,
      rule2<Graph,Tree,leaf_info<Tree>>,
      rule3<Graph,Tree,leaf_info<Tree>>,
      rule4<Graph,Tree,leaf_info<Tree>>,
      rule5<Graph,Tree,leaf_info<Tree>>,
      rule6<Graph,Tree,leaf_info<Tree>>,
    };
    int i = 0;
    bool applied = true;
    while(applied && !info.is_path()) {
      applied = false;
      for(auto rule : rules) {
        if(rule(G, T, info)) {
          ++i;
          applied = true;
          break;
        }
      }
    }
    return i;
  }
};
