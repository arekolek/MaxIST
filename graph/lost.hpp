// (C) 2014 Arek Olek

#include <functional>
#include <unordered_map>

#include <boost/functional/hash.hpp>

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
  leaf_info(Graph const * G_, Tree const & T_) : G(G_), T(T_), n(num_vertices(T)) {
    B.resize(n);
    P.resize(n*n);
    BN.resize(n*n);
    update();
  }
  bool is_path() const {
    return L.size() == 2;
  }
  std::vector<unsigned> const & leaves() const {
    return L;
  }
  std::vector<unsigned> const & leafish() const {
    return LSH;
  }
  std::vector<unsigned> const & leafish_free() const {
    return LP;
  }
  unsigned branching(unsigned l) const {
    // b(l)
    return B[l];
  }
  unsigned parent(unsigned x, unsigned l) const {
    // x->l
    return P[l*n + x];
  }
  unsigned base(unsigned x) const {
    return next(parent(x, branch(x)), x).second;
  }
  unsigned branching_neighbor(unsigned l) const {
    // b(l)->l
    return parent(branching(l), l);
  }
  unsigned branching_neighbor(unsigned l, unsigned x) const {
    // b(l)->x
    return BN[l*n + x];
  }
  unsigned branch(unsigned x) const {
    // l: x ∈ br(l)
    return BR[x];
  }
  bool on_branch(unsigned l, unsigned x) const {
    return l == x || branching(l) == x
      || (out_degree(x, T) == 2 && !on_trunk(x) && branch(x) == l);
  }
  bool on_trunk(unsigned x) const {
    return BR[x] == -1;
  }
  bool is_short(unsigned l) const {
    return parent(branching(l), l) == l;
  }
  void update() {
    L.clear();
    BR.assign(n, -1);
    LSH.clear();
    LP.clear();
    for(auto v : range(vertices(T)))
      if(out_degree(v, T) == 1) {
        L.push_back(v);
        traverse(v, T);
      }
    if(G != NULL && L.size() > 2) {
      std::vector<bool> lp(n, false);
      for(auto l : leaves())
        lp[l] = true;
      for(auto l : leaves())
        if(!is_short(l) && edge(l, branching(l), *G).second) {
          LSH.push_back(branching_neighbor(l));
          lp[l] = false;
        }
      for(auto x : range(vertices(T)))
        if(out_degree(x, T) == 2 && !on_trunk(x) && edge(branch(x), x, *G).second && x != branch(x)) {
          LSH.push_back(parent(x, branch(x)));
          lp[branch(x)] = false;
        }
      for(unsigned l = 0; l < lp.size(); ++l)
        if(lp[l])
          LP.push_back(l);
    }
  }
  void traverse(unsigned l, Tree const & T) {
    traverse(l, l, l, T);
  }
  void traverse(unsigned l, unsigned a, unsigned b, Tree const & T) {
    do {
      std::tie(a, b) = next(a, b);
      P[l*n + b] = a;
      BR[a] = l;
    } while(out_degree(b, T) == 2);
    BR[b] = l;
    if(out_degree(b, T) > 2) {
      B[l] = b;
      for(auto v : range(adjacent_vertices(b, T)))
        if(v != a)
          traverse(l, v, b, v, T);
    }
  }
  void traverse(unsigned l, unsigned blx, unsigned a, unsigned b, Tree const & T) {
    P[l*n + b] = a;
    BN[l*n + b] = blx;
    while(out_degree(b, T) == 2) {
      std::tie(a, b) = next(a, b);
      P[l*n + b] = a;
      BN[l*n + b] = blx;
    }
    if(out_degree(b, T) > 2) {
      for(auto v : range(adjacent_vertices(b, T)))
        if(v != a)
          traverse(l, blx, b, v, T);
    }
  }
  std::pair<unsigned, unsigned> next(unsigned a, unsigned b) const {
    auto it = adjacent_vertices(b, T).first;
    return std::make_pair(b, a == *it ? *(++it) : *it);
  }
private:
  Graph const* G;
  Tree const& T;
  unsigned n;
  std::vector<unsigned> L, LSH, LP;
  std::vector<int> B, BR;
  std::vector<int> P, BN;
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
  add_edge(l1, l2, T);
  remove_edge(i.branching(l1), i.branching_neighbor(l1), T);
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

template<class Graph, class Tree, class LeafInfo>
bool ruleCycleElimination(Graph& G, Tree& T, LeafInfo& info) {
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
        auto a = x;
        auto b = info.parent(x, l);
        auto bln = info.branching_neighbor(l);
        while (b != bln) {
          if (out_degree(a, T) > 2 && out_degree(b, T) > 2) {
            add_edge(l, x, T);
            remove_edge(a, b, T);
            info.update();
            return true;
          }
          if (out_degree(a, T) > 2 && supports_other(b, l)) {
            add_edge(l, x, T);
            remove_edge(a, b, T);
            info.update();
            rule1action(b, supported_other(b, l), T, info);
            return true;
          }
          if (supports_other(a, l) && out_degree(b, T) > 2) {
            add_edge(l, x, T);
            remove_edge(a, b, T);
            info.update();
            if(a != x) rule1action(a, supported_other(a, l), T, info);
            return true;
          }
          a = b;
          b = info.parent(b, l);
        }
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
      if(!extra[blx].empty() && edge(l2, blx, G).second) {
        int l, x;
        if(extra[blx].size() == 1 && extra[blx].begin()->first == l2) continue;
        std::tie(l, x) = extra[blx].begin()->first == l2 ? *(extra[blx].begin()+1) : *extra[blx].begin();
        add_edge(l, x, T);
        remove_edge(info.branching(l), blx, T);
        info.update();
        rule1action(l2, blx, T, info);
        return true;
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
        add_edge(l, x, T);
        add_edge(l, y, T);
        remove_edge(x, y, T);
        remove_edge(l, info.branching(l), T);
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
        add_edge(a, x, T);
        add_edge(b, y, T);
        remove_edge(x, y, T);
        remove_edge(info.branching(l), info.branching_neighbor(l), T);
        info.update();
        return true;
      }
    }
  }
  return false;
}

template <class Tree, class LeafInfo>
void rule7action(unsigned l1, unsigned l2, Tree& T, LeafInfo& i) {
  add_edge(i.branching_neighbor(l1), l2, T);
  remove_edge(i.branching_neighbor(l1), i.branching(l1), T);
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
        if(info.on_trunk(x) || out_degree(x, T) > 2)
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
            unsigned bn1 = info.branching_neighbor(l1);
            unsigned bn2 = info.branching_neighbor(l2);
            add_edge(l1, x, T);
            add_edge(bn1, bn2, T);
            remove_edge(info.branching(l1), bn1, T);
            remove_edge(info.branching(l2), bn2, T);
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
            unsigned bn1 = info.branching_neighbor(l1);
            unsigned bn2 = info.branching_neighbor(l2);
            add_edge(l1, x, T);
            add_edge(bn1, bn2, T);
            remove_edge(info.branching(l1), bn1, T);
            remove_edge(info.branching(l2), bn2, T);
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
  auto l = i.branch(u), x = i.base(u);
  add_edge(l, x, T);
  remove_edge(x, u, T);
  i.update();
}

template <class Graph, class Tree, class LeafInfo>
bool rule10(Graph& G, Tree& T, LeafInfo& info) {
  for(auto l : info.leaves())
    for(auto u : info.leafish())
      if(edge(u, l, T).second && info.branch(u) != l) {
        ruleA(u, T, info);
        rule1action(l, u, T, info);
        return true;
      }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule11(Graph& G, Tree& T, LeafInfo& info) {
  for(auto u : info.leafish())
    for(auto v : info.leafish())
      if(edge(u, v, T).second && info.branch(u) != info.branch(v)) {
        ruleA(u, T, info);
        ruleA(v, T, info);
        rule1action(u, v, T, info);
        return true;
      }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule12(Graph& G, Tree& T, LeafInfo& info) {
  for(auto l : info.leaves())
    for(auto u : info.leafish())
      if(edge(u, info.branching_neighbor(l), T).second
          && info.branch(u) != l) {
        ruleA(u, T, info);
        rule7action(l, u, T, info);
        return true;
      }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule13(Graph& G, Tree& T, LeafInfo& info) {
  auto const & lp = info.leafish_free();
  std::unordered_map<unsigned, unsigned> lookup;
  for(unsigned i = 0; i < lp.size(); ++i)
    lookup[lp[i]] = i;
  unsigned n = lp.size();
  std::vector<int> m(n * n, -1);
  for(unsigned i = 0; i < n; ++i)
    for(auto x : range(adjacent_vertices(lp[i], G)))
      if(!info.on_trunk(x) && info.branch(x) != lp[i])
        try {
          m[i*n + lookup.at(info.branch(x))] = x;
        } catch (std::out_of_range& e) {
        }
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < n; ++j) {
      if (m[i * n + j] >= 0 && m[j * n + i] >= 0) {
        auto l1 = lp[i];
        auto x = m[i * n + j];
        add_edge(l1, x, T);
        remove_edge(x, info.parent(x, l1), T);
        info.update();
        return true;
      }
    }
  }
  return false;
}

template <class Graph, class Tree, class LeafInfo>
bool rule14(Graph& G, Tree& T, LeafInfo& info) {
  for(auto l1 : info.leafish_free())
    for(auto l2 : info.leafish_free())
      if(l1 != l2
          && info.branching(l1) == info.branching(l2)
          && out_degree(info.branching(l1), T) == 3
          && edge(info.branching_neighbor(l1), info.branching_neighbor(l2), G).second) {
        add_edge(info.branching_neighbor(l1), info.branching_neighbor(l2), T);
        remove_edge(info.branching(l2), info.branching_neighbor(l2), T);
        info.update();
        return true;
      }
  return false;
}

template<class Graph, class Tree>
std::function<std::string(Graph&,Tree&)> make_improvement(std::string name) {
  std::vector<std::function<bool(Graph&,Tree&,leaf_info<Graph,Tree>&)>> typedef Rules;
  Rules rules = {
          rule0<Graph,Tree,leaf_info<Graph,Tree>>,
          rule1<Graph,Tree,leaf_info<Graph,Tree>>,
          rule2<Graph,Tree,leaf_info<Graph,Tree>>,
          rule3<Graph,Tree,leaf_info<Graph,Tree>>,
          rule4<Graph,Tree,leaf_info<Graph,Tree>>,
          rule5<Graph,Tree,leaf_info<Graph,Tree>>,
          ruleCycleElimination<Graph,Tree,leaf_info<Graph,Tree>>,
          rule6<Graph,Tree,leaf_info<Graph,Tree>>,
          rule6extended<Graph,Tree,leaf_info<Graph,Tree>>,
          rule7<Graph,Tree,leaf_info<Graph,Tree>>,
          rule8<Graph,Tree,leaf_info<Graph,Tree>>,
          rule9<Graph,Tree,leaf_info<Graph,Tree>>,
          rule10<Graph,Tree,leaf_info<Graph,Tree>>,
          rule11<Graph,Tree,leaf_info<Graph,Tree>>,
          rule12<Graph,Tree,leaf_info<Graph,Tree>>,
          rule13<Graph,Tree,leaf_info<Graph,Tree>>,
          rule14<Graph,Tree,leaf_info<Graph,Tree>>,
        };
  std::array<bool, 17> active;
  bool extended = false;
  if (name == "prieto")
    active = {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  else if (name == "lost-light")
    active = {0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0};
  else if (name == "lost") {
    extended = true;
    active = {0,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1};
  }
  else if (name == "lost15") {
    extended = true;
    active = {1,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1};
  }
  else if (name == "lost-ex") {
    extended = true;
    active = {0,1,0,0,0,0,1,0,1,1,1,1,1,1,1,1,1};
  }
  else if (name == "none")
    return [](Graph& G, Tree& T) { return ""; };
  else
    throw std::invalid_argument("Unknown construction method: " + name);
  return [rules, active, extended](Graph& G, Tree& T) {
    leaf_info<Graph,Tree> info(extended ? &G : NULL, T);
    bool applied = true;
    std::array<unsigned, 17> counter = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //show("step" + std::to_string(i) + ".dot", G, T);
    while(applied && !info.is_path()) {
      applied = false;
      for(unsigned i = 0; i < rules.size(); ++i) {
        if(active[i] && rules[i](G, T, info)) {
          applied = true;
          ++counter[i];
          //std::cerr << ("rule " + std::to_string(k) + "\n");
          //show("step" + std::to_string(i) + ".dot", G, T);
          break;
        }
      }
    }
    return join(counter, "-");
  };
}

