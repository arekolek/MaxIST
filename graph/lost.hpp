// (C) 2014 Arek Olek

#include <functional>
#include <unordered_map>

#include <boost/functional/hash.hpp>

#include "algorithm.hpp"
#include "debug.hpp"
#include "range.hpp"

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

template <class Graph, class Tree>
class leaf_info {
public:
  leaf_info(Graph const & g, Tree const & t, bool leafish, bool lazy)
      : _g(g),
        _t(t),
        _n(num_vertices(_t)),
        _needs_leafish(leafish),
        _lazy(lazy) {
    _branching.resize(_n);
    _parent.resize(_n * _n);
    _branching_neighbor.resize(_n * _n);
    update();
  }
  bool is_path() const {
    return _leaves.size() == 2;
  }
  std::vector<unsigned> const & leaves() const {
    return _leaves;
  }
  std::vector<unsigned> const & leafish() {
    update_leafish();
    return _leafish;
  }
  std::vector<unsigned> const & leafish_free() {
    update_leafish();
    return _leafish_free;
  }

  // b(l)
  unsigned branching(unsigned l) {
    traverse(l);
    return _branching[l];
  }
  // x->l
  unsigned parent(unsigned x, unsigned l) {
    traverse(l);
    return _parent[l*_n + x];
  }
  // b(l)->l
  unsigned branching_neighbor(unsigned l) {
    traverse(l);
    return parent(branching(l), l);
  }
  // b(l)->x
  unsigned branching_neighbor(unsigned l, unsigned x) {
    traverse(l);
    return _branching_neighbor[l*_n + x];
  }
  bool is_short(unsigned l) {
    traverse(l);
    return parent(branching(l), l) == l;
  }

  // l: x ∈ br(l)
  unsigned branch(unsigned x) {
    assert(out_degree(x, _t) < 3);
    traverse_all(); // this is not lazy, but this function almost doesn't get called
    assert(_branch[x] > -1);
    return _branch[x];
  }
  unsigned base(unsigned x) {
    assert(out_degree(x, _t) < 3);
    // traverse_all() is called in branch(x) if needed
    return next(parent(x, branch(x)), x).second;
  }
  bool on_branch(unsigned l, unsigned x) {
    // traverse(l) is called in branching(l) if needed
    return x == l
        || x == branching(l)
        || (out_degree(x, _t) == 2 && _branch[x] == (int)l); // we can access _branch directly in this case
  }
  bool on_trunk(unsigned x) {
    traverse_all();
    return _branch[x] == -1;
  }

  void update() {
    _leaves.clear();
    _traversed_all = false;
    _traversed.assign(_n, false);
    _branch.assign(_n, -1);
    _leafish.clear();
    _leafish_free.clear();
    for(auto v : range(vertices(_t)))
      if(out_degree(v, _t) == 1)
        _leaves.push_back(v);

    if(!_lazy) {
      traverse_all();
      if(_needs_leafish) update_leafish();
    }
  }

protected:
  void update_leafish() {
    if (_leafish.empty() && _leafish_free.empty()) {
      std::vector<bool> lp(_n, false);
      for (auto l : leaves())
        lp[l] = !is_short(l);
      for (auto l : leaves())
        if (!is_short(l) && edge(l, branching(l), _g).second) {
          _leafish.push_back(branching_neighbor(l));
          lp[l] = false;
        }
      for (auto x : range(vertices(_t)))
        if (out_degree(x, _t) == 2 && !on_trunk(x)
            && edge(branch(x), x, _g).second && !edge(branch(x), x, _t).second) {
          _leafish.push_back(parent(x, branch(x)));
          lp[branch(x)] = false;
        }
      for (unsigned l = 0; l < lp.size(); ++l)
        if (lp[l])
          _leafish_free.push_back(l);
    }
  }
  void traverse_all() {
    if(!_traversed_all) {
      _traversed_all = true;
      for(auto l : _leaves) traverse(l);
    }
  }
  void traverse(unsigned l) {
    if(!_traversed[l]) {
      _traversed[l] = true;
      traverse(l, l, l);
    }
  }
  void traverse(unsigned l, unsigned a, unsigned b) {
    do {
      std::tie(a, b) = next(a, b);
      _parent[l*_n + b] = a;
      _branch[a] = l;
    } while(out_degree(b, _t) == 2);
    _branch[b] = l;
    if(out_degree(b, _t) > 2) {
      _branching[l] = b;
      for(auto v : range(adjacent_vertices(b, _t)))
        if(v != a)
          traverse(l, v, b, v);
    }
  }
  void traverse(unsigned l, unsigned blx, unsigned a, unsigned b) {
    _parent[l*_n + b] = a;
    _branching_neighbor[l*_n + b] = blx;
    while(out_degree(b, _t) == 2) {
      std::tie(a, b) = next(a, b);
      _parent[l*_n + b] = a;
      _branching_neighbor[l*_n + b] = blx;
    }
    if(out_degree(b, _t) > 2) {
      for(auto v : range(adjacent_vertices(b, _t)))
        if(v != a)
          traverse(l, blx, b, v);
    }
  }
  std::pair<unsigned, unsigned> next(unsigned a, unsigned b) const {
    auto it = adjacent_vertices(b, _t).first;
    return std::make_pair(b, a == *it ? *(++it) : *it);
  }

private:
  Graph const& _g;
  Tree const& _t;
  unsigned _n;
  bool _needs_leafish;
  bool _lazy;
  bool _traversed_all;
  std::vector<bool> _traversed;
  std::vector<unsigned> _leaves, _leafish, _leafish_free;
  std::vector<int> _branching, _branch;
  std::vector<int> _parent, _branching_neighbor;
};

template <class Graph, class Tree, class LeafInfo>
bool rule0(Graph& G, Tree& T, LeafInfo& info) {
  for(auto l1 : info.leaves())
    for(auto l2 : info.leaves())
      if(edge(l1, l2, G).second) {
        auto x = l1;
        auto y = *adjacent_vertices(l1, T).first;
        while(y != l2) {
          x = y;
          y = info.parent(y, l2);
          if(out_degree(x, T) > 2 && out_degree(y, T) > 2) {
            add_edge(l1, l2, T);
            remove_edge(x, y, T);
            info.update();
            return true;
          }
        }
      }
  return false;
}

template <class Tree, class LeafInfo>
void rule1action(unsigned l1, unsigned l2, Tree& T, LeafInfo& i) {
  auto x = i.branching(l1);
  auto y = i.branching_neighbor(l1);
  add_edge(l1, l2, T);
  remove_edge(x, y, T);
  i.update();
}

template <class Graph, class Tree, class LeafInfo>
bool rule1(Graph& G, Tree& T, LeafInfo& info) {
  for(auto l1 : info.leaves())
    for(auto l2 : info.leaves())
      if(edge(l1, l2, G).second) {
        rule1action(l1, l2, T, info);
        return true;
      }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule2(Graph& G, Tree& T, LeafInfo& info) {
  for(auto l : info.leaves())
    for(auto x : range(adjacent_vertices(l, G)))
      if(!info.on_branch(l, x)) {
        auto xl = info.parent(x, l);
        if(out_degree(xl, T) > 2) {
          add_edge(l, x, T);
          remove_edge(x, xl, T);
          info.update();
          return true;
        }
      }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule3(Graph& G, Tree& T, LeafInfo& info) {
  int n = num_vertices(G);
  std::vector<std::vector<std::pair<int,int>>> extra(n);
  for(auto l : info.leaves())
    for(auto x : range(adjacent_vertices(l, G)))
      if(!info.on_branch(l, x)) {
        auto xl = info.parent(x, l);
        if(out_degree(xl, T) == 2)
          extra[xl].emplace_back(l, x);
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
        rule1action(l2, xl, T, info);
        return true;
      }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule4(Graph& G, Tree& T, LeafInfo& info) {
  for(auto l : info.leaves())
    for(auto x : range(adjacent_vertices(l, G)))
      if(!info.on_branch(l, x)) {
        auto bl = info.branching(l);
        auto blx = info.branching_neighbor(l, x);
        if(out_degree(blx, T) > 2) {
          add_edge(l, x, T);
          remove_edge(bl, blx, T);
          info.update();
          return true;
        }
      }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule5(Graph& G, Tree& T, LeafInfo& info) {
  int n = num_vertices(G);
  std::vector<std::vector<std::pair<int,int>>> extra(n);
  for(auto l : info.leaves())
    for(auto x : range(adjacent_vertices(l, G)))
      if(!info.on_branch(l, x)) {
        auto blx = info.branching_neighbor(l, x);
        if(out_degree(blx, T) == 2) {
          extra[blx].emplace_back(l, x);
        }
      }

  for(int l2 : info.leaves())
    for(int blx = 0; blx < n; ++blx)
      if(!extra[blx].empty() && edge(l2, blx, G).second && !edge(l2, blx, T).second) {
        int l1 = -1, x;
        for(auto e : extra[blx]) {
          if(e.first != l2 && e.second != l2 && e.second != blx) {
            std::tie(l1, x) = e;
            break;
          }
        }
        if(l1 == -1) continue;
        auto bl = info.branching(l1);
        add_edge(l1, x, T);
        remove_edge(bl, blx, T);
        info.update();
        rule1action(l2, blx, T, info);
        return true;
      }
  return false;
}


template<class Graph, class Tree, class LeafInfo>
bool ruleCycleElimination1(Graph& G, Tree& T, LeafInfo& info) {
  for (auto l : info.leaves())
    for (auto x : range(adjacent_vertices(l, G)))
      if (!info.on_branch(l, x)) {
        auto a = info.parent(x, l);
        auto b = info.parent(a, l);
        auto bl = info.branching(l);
        if(a == bl) continue;
        while (b != bl) {
          if (out_degree(a, T) > 2 && out_degree(b, T) > 2) {
            add_edge(l, x, T);
            remove_edge(a, b, T);
            info.update();
            return true;
          }
          a = b;
          b = info.parent(b, l);
        }
      }
  return false;
}


template<class Graph, class Tree, class LeafInfo>
bool ruleCycleElimination2(Graph& G, Tree& T, LeafInfo& info) {
  int n = num_vertices(G);
  std::vector<std::vector<uint>> supported(n);
  for (auto l : info.leaves())
    for (auto x : range(adjacent_vertices(l, G)))
      if (!info.on_branch(l, x))
        supported[x].push_back(l);

  auto supports_other = [supported](uint x, uint l){
    return supported[x].size() > 1 || (supported[x].size() == 1 && *supported[x].begin() != l);
  };

  auto supported_other = [supported](uint x, uint l){
    auto it = supported[x].begin();
    return *it == l ? *(it+1) : *it;
  };

  for (auto l : info.leaves())
    for (auto x : range(adjacent_vertices(l, G)))
      if (!info.on_branch(l, x)) {
        auto a = info.parent(x, l);
        auto b = info.parent(a, l);
        auto bl = info.branching(l);
        if(a == bl) continue;
        while (b != bl) {
          if (supports_other(a, l) && out_degree(a, T) == 2 && out_degree(b, T) > 2) std::swap(a, b);
          if (out_degree(a, T) > 2 && out_degree(b, T) == 2 && supports_other(b, l)) {
            add_edge(l, x, T);
            remove_edge(a, b, T);
            info.update();
            rule1action(b, supported_other(b, l), T, info);
            return true;
          }
          a = b;
          b = info.parent(b, l);
        }
      }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule6(Graph& G, Tree& T, LeafInfo& info) {
  for(auto e : range(edges(T))) {
    auto x = source(e, T);
    auto y = target(e, T);
    for(auto l : info.leaves()) {
      if(info.is_short(l)
        && !edge(l, x, T).second && !edge(l, y, T).second
        && edge(l, x, G).second && edge(l, y, G).second) {
        auto bl = info.branching(l);
        add_edge(l, x, T);
        add_edge(l, y, T);
        remove_edge(x, y, T);
        remove_edge(l, bl, T);
        info.update();
        return true;
      }
    }
  }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule6extended(Graph& G, Tree& T, LeafInfo& info) {
  for(auto e : range(edges(T))) {
    auto x = source(e, T);
    auto y = target(e, T);
    auto check = [&](uint a, uint b) {
      return !edge(a, x, T).second && !edge(b, y, T).second
      && edge(a, x, G).second && edge(b, y, G).second;
    };
    for(auto l : info.leaves()) {
      if(info.on_branch(l, x) && x != info.branching(l)) continue;
      if(info.on_branch(l, y) && y != info.branching(l)) continue;
      auto a = l;
      auto b = info.branching_neighbor(l);
      if(check(a, b) || check(b, a)) {
        if(check(b, a)) std::swap(a, b);
        auto bl = info.branching(l);
        auto bln = info.branching_neighbor(l);
        add_edge(a, x, T);
        add_edge(b, y, T);
        remove_edge(x, y, T);
        remove_edge(bl, bln, T);
        info.update();
        return true;
      }
    }
  }
  return false;
}

template <class Tree, class LeafInfo>
void rule7action(unsigned l1, unsigned l2, Tree& T, LeafInfo& i) {
  auto bl = i.branching(l1);
  auto bln = i.branching_neighbor(l1);
  add_edge(bln, l2, T);
  remove_edge(bln, bl, T);
  i.update();
}

template <class Graph, class Tree, class LeafInfo>
bool rule7(Graph& G, Tree& T, LeafInfo& info) {
  for(auto l1 : info.leaves()) if(!info.is_short(l1)) {
    for(auto l2 : info.leaves()) if(/*!info.is_short(l2) && */l1 != l2) {
      if(edge(info.branching_neighbor(l1), l2, G).second) {
        rule7action(l1, l2, T, info);
        return true;
      }
    }
  }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule8(Graph& G, Tree& T, LeafInfo& info) {
  std::vector<unsigned> lg;
  for(auto l : info.leaves())
    if(!info.is_short(l))
      lg.push_back(l);
  unsigned n = lg.size();
  std::vector<bool> m(n * n);
  for (unsigned i = 0; i < n; ++i) {
    unsigned l1 = lg[i];
    for (unsigned j = 0; j < n; ++j) {
      unsigned l2 = lg[j];
      m[i*n + j] =
        info.branching(l1) != info.branching(l2) &&
        edge(
          info.branching_neighbor(l1),
          info.branching_neighbor(l2),
          G).second;
    }
  }
  for (unsigned i = 0; i < n; ++i) {
    unsigned l1 = lg[i];
    unsigned count = 0, l2 = 0;
    bool ok = false;
    for(auto x : range(adjacent_vertices(l1, G)))
      if(!info.on_branch(l1, x)) {
        ++count;
        if(info.on_trunk(x) || out_degree(x, T) > 2) // FIXME why "|| out_degree(x, T) > 2"
          ok = true;
        else if(count == 1)
          l2 = info.branch(x);
        else if(l2 != info.branch(x))
          ok = true;
        if(ok) break;
      }
    if(count == 0)
      for(unsigned j = 0; j < n; ++j)
        m[i*n + j] = false;
    else if(!ok)
      for(unsigned j = 0; j < n; ++j)
        if(lg[j] == l2)
          m[i*n + j] = false;
  }
  for(unsigned i = 0; i < n; ++i) {
    unsigned l1 = lg[i];
    for(unsigned j = 0; j < n; ++j) {
      unsigned l2 = lg[j];
      if(m[i*n + j]) {
        for(auto x : range(adjacent_vertices(l1, G)))
          if(!info.on_branch(l1, x) && (!info.on_branch(l2, x) || out_degree(x, T) > 2)) {
            auto b1 = info.branching(l1);
            auto b2 = info.branching(l2);
            auto bn1 = info.branching_neighbor(l1);
            auto bn2 = info.branching_neighbor(l2);
            add_edge(l1, x, T);
            add_edge(bn1, bn2, T);
            remove_edge(b1, bn1, T);
            remove_edge(b2, bn2, T);
            info.update();
            return true;
          }
        assert(false);
      }
    }
  }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule9(Graph& G, Tree& T, LeafInfo& info) {
  std::vector<unsigned> lg;
  for(auto l : info.leaves())
    if(!info.is_short(l))
      lg.push_back(l);
  unsigned n = lg.size();
  std::vector<bool> m(n * n);
  for (unsigned i = 0; i < n; ++i) {
    unsigned l1 = lg[i];
    for (unsigned j = 0; j < n; ++j) {
      unsigned l2 = lg[j];
      m[i*n + j] =
        info.branching(l1) == info.branching(l2) &&
        out_degree(info.branching(l1), T) >= 4 &&
        edge(
          info.branching_neighbor(l1),
          info.branching_neighbor(l2),
          G).second;
    }
  }
  for (unsigned i = 0; i < n; ++i) {
    unsigned l1 = lg[i];
    unsigned count = 0, l2 = 0;
    bool ok = false;
    for(auto x : range(adjacent_vertices(l1, G)))
      if(!info.on_branch(l1, x)) {
        if(info.on_trunk(x) || (out_degree(x, T) > 2 && x != info.branching(l1)))
          ok = true;
        else if(out_degree(x, T) == 2) {
          ++count;
          if(count == 1)
            l2 = info.branch(x);
          else if(l2 != info.branch(x))
            ok = true;
        }
        if(ok) break;
      }
    if(!ok)
      for(unsigned j = 0; j < lg.size(); ++j)
        if(count == 0 || lg[j] == l2)
          m[i*n + j] = false;
  }
  for(unsigned i = 0; i < lg.size(); ++i) {
    unsigned l1 = lg[i];
    for(unsigned j = 0; j < lg.size(); ++j) {
      unsigned l2 = lg[j];
      if(m[i*n + j]) {
        for(auto x : range(adjacent_vertices(l1, G)))
          if(!info.on_branch(l1, x) && !info.on_branch(l2, x)) {
            auto b1 = info.branching(l1);
            auto b2 = info.branching(l2);
            auto bn1 = info.branching_neighbor(l1);
            auto bn2 = info.branching_neighbor(l2);
            add_edge(l1, x, T);
            add_edge(bn1, bn2, T);
            remove_edge(b1, bn1, T);
            remove_edge(b2, bn2, T);
            info.update();
            return true;
          }
        assert(false);
      }
    }
  }
  return false;
}

template <class Tree, class LeafInfo>
void ruleA(unsigned u, Tree& T, LeafInfo& i) {
  auto l = i.branch(u);
  auto x = i.base(u);
  add_edge(l, x, T);
  remove_edge(x, u, T);
  i.update();
}

template <class Graph, class Tree, class LeafInfo>
bool rule10(Graph& G, Tree& T, LeafInfo& info) {
  for(auto l2 : info.leaves())
    for(auto u : info.leafish())
      if(edge(u, l2, G).second && info.branch(u) != l2) {
        ruleA(u, T, info);
        rule1action(l2, u, T, info);
        return true;
      }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule11(Graph& G, Tree& T, LeafInfo& info) {
  for(auto u : info.leafish())
    for(auto v : info.leafish())
      if(edge(u, v, G).second && info.branch(u) != info.branch(v)) {
        ruleA(u, T, info);
        ruleA(v, T, info);
        rule1action(u, v, T, info);
        return true;
      }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule12(Graph& G, Tree& T, LeafInfo& info) {
  for(auto l2 : info.leaves())
    for(auto u : info.leafish())
      if(edge(u, info.branching_neighbor(l2), G).second
          && info.branch(u) != l2) {
        ruleA(u, T, info);
        rule7action(l2, u, T, info);
        return true;
      }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule13(Graph& G, Tree& T, LeafInfo& info) {
  auto const & lp = info.leafish_free();
  unsigned n = lp.size();
  std::vector<int> m(n * n, -1);
  for(unsigned i = 0; i < n; ++i)
    for(unsigned j = 0; j < n; ++j)
      if(i != j)
        for(auto x : range(adjacent_vertices(lp[i], G)))
          if(out_degree(x, T) == 2 && info.on_branch(lp[j], x))
              m[i*n + j] = x;
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < n; ++j) {
      if (m[i * n + j] >= 0 && m[j * n + i] >= 0) {
        auto l1 = lp[i];
        auto x = m[i * n + j];
        auto xl = info.parent(x, l1);
        add_edge(l1, x, T);
        remove_edge(x, xl, T);
        info.update();
        return true;
      }
    }
  }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule14(Graph& G, Tree& T, LeafInfo& info) {
  if(info.leaves().size() > 3)
    for(auto l1 : info.leafish_free())
      for(auto l2 : info.leafish_free())
        if(l1 != l2
            && info.branching(l1) == info.branching(l2)
            && out_degree(info.branching(l1), T) == 3
            && edge(info.branching_neighbor(l1), info.branching_neighbor(l2), G).second) {
          auto bn1 = info.branching_neighbor(l1);
          auto bn2 = info.branching_neighbor(l2);
          auto b2 = info.branching(l2);
          add_edge(bn1, bn2, T);
          remove_edge(b2, bn2, T);
          info.update();
          return true;
        }
  return false;
}

template<class Graph, class Tree>
std::function<std::vector<unsigned>(Graph&,Tree&)> make_improvement(std::string name) {
  std::vector<std::function<bool(Graph&,Tree&,leaf_info<Graph,Tree>&)>> typedef Rules;
  Rules rules = {
          rule0<Graph,Tree,leaf_info<Graph,Tree>>,
          rule1<Graph,Tree,leaf_info<Graph,Tree>>,
          rule2<Graph,Tree,leaf_info<Graph,Tree>>,
          rule3<Graph,Tree,leaf_info<Graph,Tree>>,
          rule4<Graph,Tree,leaf_info<Graph,Tree>>,
          rule5<Graph,Tree,leaf_info<Graph,Tree>>,
          rule6<Graph,Tree,leaf_info<Graph,Tree>>,
          rule7<Graph,Tree,leaf_info<Graph,Tree>>,
          rule8<Graph,Tree,leaf_info<Graph,Tree>>,
          rule9<Graph,Tree,leaf_info<Graph,Tree>>,
          rule10<Graph,Tree,leaf_info<Graph,Tree>>,
          rule11<Graph,Tree,leaf_info<Graph,Tree>>,
          rule12<Graph,Tree,leaf_info<Graph,Tree>>,
          rule13<Graph,Tree,leaf_info<Graph,Tree>>,
          rule14<Graph,Tree,leaf_info<Graph,Tree>>,
          ruleCycleElimination1<Graph,Tree,leaf_info<Graph,Tree>>,
          ruleCycleElimination2<Graph,Tree,leaf_info<Graph,Tree>>,
          rule6extended<Graph,Tree,leaf_info<Graph,Tree>>,
        };

  std::vector<unsigned> active;
  bool lazy = name.find("lazy") != std::string::npos;
  bool rand = name.find("rand") != std::string::npos;
  name = name.substr(0, name.find("/"));

  if (name == "prieto") {
    active = {1};
  }
  else if (name == "lost-light") {
    active = {1,2,3,4,5};
  }
  else if (name == "lost") {
    active = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
  }
  else if (name == "lost15") {
    active = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
  }
  else if (name == "lost-ex") {
    active = {1,2,3,4,5,15,16,6,17,7,8,9,10,11,12,13,14};
  }
  else if (name == "none") {
    return [active](Graph& G, Tree& T) {
      return std::vector<unsigned>(1, 0);
    };
  }
  else {
    std::stringstream ss(name);
    unsigned i;
    while (ss >> i) {
      active.push_back(i);
      if (ss.peek() == '+') ss.ignore();
    }
  }

  bool leafish = false;
  for(auto rule : active)
    if(9 > rule && rule < 15) leafish = true;

  return [=](Graph& G, Tree& T) mutable {
    leaf_info<Graph,Tree> info(G, T, leafish, lazy);
    bool applied = true;
    std::vector<unsigned> counter(active.size(), 0);
    auto order = active;
#ifndef NDEBUG
    unsigned steps = 0;
    //show("tree-" + std::to_string(steps) + ".dot", G, T);
#endif
    while(applied && !info.is_path()) {
      assert(++steps < num_vertices(G));
      assert(num_edges(T) == num_vertices(T)-1);
      assert(is_connected(T));
      applied = false;
      if(rand) std::random_shuffle(order.begin(), order.end());
      for(unsigned i : order) {
        if(rules[i](G, T, info)) {
          applied = true;
          ++counter[find_index(active, i)];
#ifndef NDEBUG
          //show("tree-" + std::to_string(steps) + "-" + std::to_string(i) + ".dot", G, T);
#endif
          break;
        }
      }
    }
    return counter;
  };
}

