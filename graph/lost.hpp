// (C) 2014 Arek Olek

#include <functional>
#include <tuple>

#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/functional/hash.hpp>
#include <boost/optional.hpp>
#include <boost/range/algorithm/find_if.hpp>

#include "algorithm.hpp"
#include "debug.hpp"
#include "range.hpp"

using std::vector;

using boost::find_if;
using boost::optional;

typedef std::pair<unsigned, unsigned> Edge;


template<class Graph, class Tree>
class leaf_info {
 public:
  leaf_info(Graph const & g, Tree const & t, bool leafish, bool lazy)
      : graph_(g),
        tree_(t),
        num_vertices_(num_vertices(tree_)),
        needs_leafish_(leafish),
        lazy_(lazy) {
    branching_.resize(num_vertices_);
    parent_.resize(num_vertices_ * num_vertices_);
    branching_neighbor_.resize(num_vertices_ * num_vertices_);
    update();
  }
  bool is_path() const {
    return leaves_.size() == 2;
  }
  vector<unsigned> const & leaves() const {
    return leaves_;
  }
  vector<unsigned> const & leafish() {
    update_leafish();
    return leafish_;
  }
  vector<unsigned> const & leafish_free() {
    update_leafish();
    return leafish_free_;
  }

  // b(l)
  unsigned branching(unsigned l) {
    traverse(l);
    return branching_[l];
  }
  // x->l
  unsigned parent(unsigned x, unsigned l) {
    traverse(l);
    return parent_[l * num_vertices_ + x];
  }
  // b(l)->l
  unsigned branching_neighbor(unsigned l) {
    traverse(l);
    return parent(branching(l), l);
  }
  // b(l)->x
  unsigned branching_neighbor(unsigned l, unsigned x) {
    traverse(l);
    return branching_neighbor_[l * num_vertices_ + x];
  }
  bool is_short(unsigned l) {
    traverse(l);
    return parent(branching(l), l) == l;
  }
  bool is_long(unsigned l) {
    traverse(l);
    return parent(branching(l), l) != l;
  }

  // l: x ∈ br(l)
  unsigned branch(unsigned x) {
    assert(out_degree(x, tree_) < 3);
    traverse_all();  // this is not lazy, but this function almost doesn't get called
    assert(branch_[x] > -1);
    return branch_[x];
  }
  unsigned base(unsigned x) {
    assert(out_degree(x, tree_) < 3);
    // traverse_all() is called in branch(x) if needed
    return next(parent(x, branch(x)), x).second;
  }
  bool on_branch(unsigned l, unsigned x) {
    // traverse(l) is called in branching(l) if needed
    return x == l
        || x == branching(l)
        || (out_degree(x, tree_) == 2 && branch_[x] == (int) l);
  }
  bool on_trunk(unsigned x) {
    traverse_all();
    return branch_[x] == -1;
  }

  auto support(unsigned l) {
    auto neighbors = adjacent_vertices(l, graph_);
    auto supports = [=](unsigned x) {return !on_branch(l, x);};
    auto start = filter(supports, neighbors.first, neighbors.second);
    auto end = filter(supports, neighbors.second, neighbors.second);
    return range(start, end);
  }

  void update() {
    leaves_.clear();
    traversed_all_ = false;
    traversed_.assign(num_vertices_, false);
    branch_.assign(num_vertices_, -1);
    leafish_.clear();
    leafish_free_.clear();
    for (auto v : range(vertices(tree_)))
      if (out_degree(v, tree_) == 1) leaves_.push_back(v);

    if (!lazy_) {
      traverse_all();
      if (needs_leafish_) update_leafish();
    }
  }

 protected:
  void update_leafish() {
    if (leafish_.empty() && leafish_free_.empty()) {
      vector<bool> lp(num_vertices_, false);
      for (auto l : leaves()) lp[l] = is_long(l);
      for (auto l : leaves())
        if (is_long(l) && edge(l, branching(l), graph_).second) {
          leafish_.push_back(branching_neighbor(l));
          lp[l] = false;
        }
      for (auto x : range(vertices(tree_)))
        if (out_degree(x, tree_) == 2 && !on_trunk(x)
            && edge(branch(x), x, graph_).second
            && !edge(branch(x), x, tree_).second) {
          leafish_.push_back(parent(x, branch(x)));
          lp[branch(x)] = false;
        }
      for (unsigned l = 0; l < lp.size(); ++l)
        if (lp[l]) leafish_free_.push_back(l);
    }
  }
  void traverse_all() {
    if (!traversed_all_) {
      traversed_all_ = true;
      for (auto l : leaves_) traverse(l);
    }
  }
  void traverse(unsigned l) {
    if (!traversed_[l]) {
      traversed_[l] = true;
      traverse(l, l, l);
    }
  }
  void traverse(unsigned l, unsigned a, unsigned b) {
    do {
      std::tie(a, b) = next(a, b);
      parent_[l * num_vertices_ + b] = a;
      branch_[a] = l;
    } while (out_degree(b, tree_) == 2);
    branch_[b] = l;
    if (out_degree(b, tree_) > 2) {
      branching_[l] = b;
      for (auto v : range(adjacent_vertices(b, tree_)))
        if (v != a) traverse(l, v, b, v);
    }
  }
  void traverse(unsigned l, unsigned blx, unsigned a, unsigned b) {
    parent_[l * num_vertices_ + b] = a;
    branching_neighbor_[l * num_vertices_ + b] = blx;
    while (out_degree(b, tree_) == 2) {
      std::tie(a, b) = next(a, b);
      parent_[l * num_vertices_ + b] = a;
      branching_neighbor_[l * num_vertices_ + b] = blx;
    }
    if (out_degree(b, tree_) > 2) {
      for (auto v : range(adjacent_vertices(b, tree_)))
        if (v != a) traverse(l, blx, b, v);
    }
  }
  Edge next(unsigned a, unsigned b) const {
    auto it = adjacent_vertices(b, tree_).first;
    return std::make_pair(b, a == *it ? *(++it) : *it);
  }

 private:
  Graph const& graph_;
  Tree const& tree_;
  unsigned num_vertices_;
  bool needs_leafish_;
  bool lazy_;
  bool traversed_all_;
  vector<bool> traversed_;
  vector<unsigned> leaves_, leafish_, leafish_free_;
  vector<int> branching_, branch_;
  vector<int> parent_, branching_neighbor_;
};

template<class Graph, class Tree, class LeafInfo>
bool rule0(Graph& G, Tree& T, LeafInfo& info) {
  auto next = [&](uint x, uint l) {return info.parent(x, l);};
  for (auto l1 : info.leaves())
    for (auto l2 : info.leaves())
      if (edge(l1, l2, G).second)
        for (auto x = l1, y = next(x, l2); y != l2; x = y, y = next(y, l2))
          if (out_degree(x, T) > 2 && out_degree(y, T) > 2) {
            add_edge(l1, l2, T);
            remove_edge(x, y, T);
            info.update();
            return true;
          }
  return false;
}

template<class Tree, class LeafInfo>
void rule1action(unsigned l1, unsigned l2, Tree& T, LeafInfo& i) {
  assert(out_degree(l1, T) == 1);
  assert(out_degree(l2, T) == 1);
  if (i.is_path()) return;
  auto x = i.branching(l1);
  auto y = i.branching_neighbor(l1);
  add_edge(l1, l2, T);
  remove_edge(x, y, T);
  i.update();
}

// Eliminates a pair of conflicted leaves, attaching one branch to the other.
template<class Graph, class Tree, class LeafInfo>
bool rule1(Graph& G, Tree& T, LeafInfo& info) {
  for (auto l1 : info.leaves())
    for (auto l2 : info.leaves())
      if (edge(l1, l2, G).second) {
        rule1action(l1, l2, T, info);
        return true;
      }
  return false;
}

// Connects a leaf l to a supporting vertex x that neighbors a branching xl.
// Removes edge (x, xl) to delete the cycle.
template<class Graph, class Tree, class LeafInfo>
bool rule2(Graph& G, Tree& T, LeafInfo& info) {
  for (auto l : info.leaves())
    for (auto x : info.support(l)) {
      auto xl = info.parent(x, l);
      if (out_degree(xl, T) > 2) {
        add_edge(l, x, T);
        remove_edge(x, xl, T);
        info.update();
        return true;
      }
    }
  return false;
}

// Connects a leaf l1 to a supporting vertex x if its forwarding neighbor xl
// has a nontree edge to another leaf l2. Removes edge (x, xl) thus creating
// a pair of conflicted leaves l2, xl, to which rule1 is then applied.
template<class Graph, class Tree, class LeafInfo>
bool rule3(Graph& G, Tree& T, LeafInfo& info) {
  vector<optional<Edge>> support_edge(num_vertices(G));
  // Find candidate triples of vertices.
  for (auto l1 : info.leaves()) {
    for (auto x : info.support(l1)) {
      auto xl = info.parent(x, l1);
      // Salamon didn't consider that l1 may be the only leaf neighbor of xl,
      // in such case the rule can't be executed, so it's necessary
      // to record which leaf is supported by x to compare it with l2 later.
      // To save on computation record x also.
      if (out_degree(xl, T) == 2) support_edge[xl] = Edge(l1, x);
    }
  }

  // Find second non-tree edge.
  for (auto l2 : info.leaves()) {
    // When using adjacency list, it's easier to just get neighbors of l2.
    // With an adjacency matrix it might be better to check adjacency
    // for each pair of leaf and candidate vertex xl.
    for (auto xl : range(adjacent_vertices(l2, G))) {
      // Checking !edge(l2, xl, T).second is not necessary:
      // - if l2 != x, then xl is a branching, so there is no support_edge[xl],
      // - if l2 == x, rule1 is applicable. This one will have the same effect.
      auto e = support_edge[xl];
      if (e && e->first != l2) {
        auto l1 = e->first, x = e->second;
        add_edge(l1, x, T);
        remove_edge(x, xl, T);
        info.update();
        rule1action(l2, xl, T, info);
        return true;
      }
    }
  }
  return false;
}

template<class Graph, class Tree, class LeafInfo>
bool rule20(Graph& G, Tree& T, LeafInfo& info) {
  // AKA naive rule 3
  uint xl;
  for (auto l1 : info.leaves()) {
    for (auto x : info.support(l1)) {
      if (out_degree(xl = info.parent(x, l1), T) == 2) {
        for (auto l2 : range(adjacent_vertices(xl, G))) {
          // No need to check x != l2, it's implied by (l2, xl) being non-tree.
          if (l1 != l2 && out_degree(l2, T) == 1 && !edge(l2, xl, T).second) {
            add_edge(l1, x, T);
            remove_edge(x, xl, T);
            info.update();
            rule1action(l2, xl, T, info);
            return true;
          }
        }
      }
    }
  }
  return false;
}

template<class Graph, class Tree, class LeafInfo>
bool rule4(Graph& G, Tree& T, LeafInfo& info) {
  for (auto l : info.leaves())
    for (auto x : info.support(l)) {
      auto bl = info.branching(l);
      auto blx = info.branching_neighbor(l, x);
      if (out_degree(blx, T) > 2) {
        add_edge(l, x, T);
        remove_edge(bl, blx, T);
        info.update();
        return true;
      }
    }
  return false;
}

template<class Graph, class Tree, class LeafInfo>
bool rule5(Graph& G, Tree& T, LeafInfo& info) {
  vector<Edge> extra[num_vertices(G)];
  for (auto l : info.leaves())
    for (auto x : info.support(l)) {
      auto blx = info.branching_neighbor(l, x);
      if (x != blx && out_degree(blx, T) == 2) extra[blx].emplace_back(l, x);
    }

  for (auto l2 : info.leaves())
    for (auto blx : range(adjacent_vertices(l2, G)))
      if (!edge(l2, blx, T).second) {
        auto e = find_if(extra[blx], [=](Edge e) {
          return e.first != l2 && e.second != l2;
        });
        assert(std::distance(extra[blx].begin(), e) == 0);
        if (e == extra[blx].end()) continue;
        auto l1 = e->first, x = e->second;
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
bool rule21(Graph& G, Tree& T, LeafInfo& info) {
  // AKA naive rule 5
  for (auto l1 : info.leaves())
    for (auto x : info.support(l1)) {
      auto blx = info.branching_neighbor(l1, x);
      if (x != blx && out_degree(blx, T) == 2)
        for (auto l2 : range(adjacent_vertices(blx, G)))
          if (l1 != l2 && x != l2 && out_degree(l2, T) == 1 && !edge(l2, blx, T).second) {
            auto bl = info.branching(l1);
            add_edge(l1, x, T);
            remove_edge(bl, blx, T);
            info.update();
            rule1action(l2, blx, T, info);
            return true;
          }
    }
  return false;
}

template<class Graph, class Tree, class LeafInfo>
bool rule6(Graph& G, Tree& T, LeafInfo& info) {
  for (auto e : range(edges(T))) {
    auto x = source(e, T), y = target(e, T);
    for (auto l : info.leaves()) {
      if (info.is_short(l)
          && !edge(l, x, T).second && !edge(l, y, T).second
          &&  edge(l, x, G).second &&  edge(l, y, G).second) {
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

template<class Tree, class LeafInfo>
void rule7action(unsigned l1, unsigned l2, Tree& T, LeafInfo& i) {
  auto bl = i.branching(l1);
  auto bln = i.branching_neighbor(l1);
  add_edge(bln, l2, T);
  remove_edge(bln, bl, T);
  i.update();
}

template<class Graph, class Tree, class LeafInfo>
bool rule7(Graph& G, Tree& T, LeafInfo& info) {
  for (auto l1 : info.leaves())
    for (auto l2 : info.leaves())
      if (l1 != l2 && edge(info.branching_neighbor(l1), l2, G).second) {
        rule7action(l1, l2, T, info);
        return true;
      }
  return false;
}

template<class Graph, class Tree, class LeafInfo, class Condition>
bool double_relocation(Graph& G, Tree& T, LeafInfo& info, Condition condition) {
  vector<unsigned> lg;
  for (auto l : info.leaves()) if (info.is_long(l)) lg.push_back(l);
  unsigned n = lg.size();
  vector<bool> m(n * n, false);
  for (unsigned i = 0; i < n * n; ++i) m[i] = condition(lg[i / n], lg[i % n]);

  for (unsigned i = 0; i < n; ++i) {
    auto l1 = lg[i];
    auto xs = info.support(l1);
    if (boost::empty(xs)) {
      for (unsigned j = 0; j < n; ++j) m[i * n + j] = false;
    } else {
      auto outside =
          [&](unsigned x) {return out_degree(x, T) > 2 || info.on_trunk(x);};
      if (boost::algorithm::any_of(xs, outside)) continue;
      auto l2 = info.branch(*xs.begin());
      if (info.is_short(l2)) continue;
      auto on_one_branch = [&](unsigned x) {return info.on_branch(l2, x);};
      if (boost::algorithm::all_of(xs, on_one_branch))
        m[i * n + find_index(lg, l2)] = false;
    }
  }

  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < n; ++j)
      if (m[i * n + j]) {
        unsigned l1 = lg[i], l2 = lg[j];
        for (auto x : info.support(l1))
          if (out_degree(x, T) > 2 || !info.on_branch(l2, x)) {
            auto b1 = info.branching(l1), bn1 = info.branching_neighbor(l1);
            auto b2 = info.branching(l2), bn2 = info.branching_neighbor(l2);
            add_edge(l1, x, T);
            add_edge(bn1, bn2, T);
            remove_edge(b1, bn1, T);
            remove_edge(b2, bn2, T);
            info.update();
            return true;
          }
        assert(false);
      }
  return false;
}

template<class Graph, class Tree, class LeafInfo>
bool rule8(Graph& G, Tree& T, LeafInfo& info) {
  return double_relocation(G, T, info, [&](auto l1, auto l2) {
    return info.branching(l1) != info.branching(l2)
        && edge(info.branching_neighbor(l1),
                info.branching_neighbor(l2),
                G).second;
  });
}

template<class Graph, class Tree, class LeafInfo>
bool rule9(Graph& G, Tree& T, LeafInfo& info) {
  return double_relocation(G, T, info, [&](auto l1, auto l2) {
    return info.branching(l1) == info.branching(l2)
        && out_degree(info.branching(l1), T) >= 4
        && edge(info.branching_neighbor(l1),
                info.branching_neighbor(l2),
                G).second;
  });
}

template<class Tree, class LeafInfo>
void ruleA(unsigned u, Tree& T, LeafInfo& i) {
  auto l = i.branch(u);
  auto x = i.base(u);
  add_edge(l, x, T);
  remove_edge(x, u, T);
  i.update();
}

template<class Graph, class Tree, class LeafInfo>
bool rule10(Graph& G, Tree& T, LeafInfo& info) {
  for (auto l : info.leaves())
    for (auto u : info.leafish())
      if (edge(u, l, G).second && !info.on_branch(l, u)) {
        ruleA(u, T, info);
        rule1action(u, l, T, info);
        return true;
      }
  return false;
}

template<class Graph, class Tree, class LeafInfo>
bool rule11(Graph& G, Tree& T, LeafInfo& info) {
  for (auto u : info.leafish())
    for (auto v : info.leafish())
      if (edge(u, v, G).second && info.branch(u) != info.branch(v)) {
        ruleA(u, T, info);
        ruleA(v, T, info);
        rule1action(u, v, T, info);
        return true;
      }
  return false;
}

template<class Graph, class Tree, class LeafInfo>
bool rule12(Graph& G, Tree& T, LeafInfo& i) {
  for (auto l : i.leaves())
    for (auto u : i.leafish())
      if (edge(u, i.branching_neighbor(l), G).second && !i.on_branch(l, u)) {
        ruleA(u, T, i);
        rule7action(l, u, T, i);
        return true;
      }
  return false;
}

template<class Graph, class Tree, class LeafInfo>
bool rule13(Graph& G, Tree& T, LeafInfo& info) {
  auto const & lp = info.leafish_free();
  unsigned n = lp.size();
  vector<int> m(n * n, -1);
  for (unsigned i = 0; i < n; ++i)
    for (auto x : info.support(lp[i])) {
      if (out_degree(x, T) != 2 || info.on_trunk(x)) continue;
      auto j = find_index(lp, info.branch(x));
      if (j < n) m[i * n + j] = x;
    }
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < n; ++j)
      if (m[i * n + j] >= 0 && m[j * n + i] >= 0) {
        auto l1 = lp[i];
        auto x = m[i * n + j];
        auto xl = info.parent(x, l1);
        add_edge(l1, x, T);
        remove_edge(x, xl, T);
        info.update();
        return true;
      }
  return false;
}

template<class Graph, class Tree, class LeafInfo>
bool rule14(Graph& G, Tree& T, LeafInfo& i) {
  if (i.leaves().size() < 4) return false;
  for (auto l1 : i.leafish_free())
    for (auto l2 : i.leafish_free()) {
      auto b = i.branching(l1);
      auto bn1 = i.branching_neighbor(l1), bn2 = i.branching_neighbor(l2);
      if (l1 != l2 && b == i.branching(l2) && out_degree(b, T) == 3
          && edge(bn1, bn2, G).second) {
        add_edge(bn1, bn2, T);
        remove_edge(b, bn2, T);
        i.update();
        return true;
      }
    }
  return false;
}

template<class Graph, class Tree, class LeafInfo>
bool rule15(Graph& G, Tree& T, LeafInfo& info) {
  for (auto l : info.leaves()) {
    auto next = [&](uint x) {return info.parent(x, l);};
    auto bl = info.branching(l);
    for (auto x : info.support(l))
      for (auto a = x, b = next(a); a != bl; a = b, b = next(b))
        if ((a == x || out_degree(a, T) > 2) && out_degree(b, T) > 2) {
          add_edge(l, x, T);
          remove_edge(a, b, T);
          info.update();
          return true;
        }
  }
  return false;
}

template<class Graph, class Tree, class LeafInfo>
bool rule16(Graph& G, Tree& T, LeafInfo& info) {
  auto is_branching = [&](uint x) {return out_degree(x, T) > 2;};
  auto supported_other = [&](uint b, uint l, uint x) -> boost::optional<uint> {
    for(uint v : range(adjacent_vertices(b, G)))
      if(v != l && v != x && out_degree(v, T) == 1)
        return v;
    return boost::none;
  };

  for (auto l : info.leaves()) {
    auto next = [&](uint x) {return info.parent(x, l);};
    auto bl = info.branching(l);
    for (auto x : info.support(l))
      for (auto a = x, b = next(a), c = next(b); b != bl; a = b, b = c, c = next(c))
        if (!is_branching(b) && (x == a || is_branching(a) || is_branching(c)))
          if (auto l2 = supported_other(b, l, x)) {
            add_edge(l, x, T);
            remove_edge(x == a || is_branching(a) ? a : c, b, T);
            info.update();
            rule1action(b, *l2, T, info);
            return true;
          }
  }
  return false;
}

template<class Graph, class Tree, class LeafInfo>
bool rule17(Graph& G, Tree& T, LeafInfo& info) {
  for (auto e : range(edges(T))) {
    auto x = source(e, T);
    auto y = target(e, T);
    auto check = [&](uint a, uint b) {
      return !edge(a, x, T).second && !edge(b, y, T).second
          &&  edge(a, x, G).second &&  edge(b, y, G).second;
    };
    for (auto l : info.leaves()) {
      if (info.on_branch(l, x) && x != info.branching(l)) continue;
      if (info.on_branch(l, y) && y != info.branching(l)) continue;
      auto a = l;
      auto b = info.branching_neighbor(l);
      if (check(a, b) || check(b, a)) {
        if (check(b, a)) std::swap(a, b);
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

template<class Graph, class Tree, class LeafInfo>
bool rule18(Graph& G, Tree& T, LeafInfo& info) {
  // AKA rules 15 and 16 in one
  auto is_branching = [&](uint x) {return out_degree(x, T) > 2;};
  auto supported_other = [&](uint b, uint l, uint x) -> boost::optional<uint> {
    for(uint v : range(adjacent_vertices(b, G)))
      if(v != l && v != x && out_degree(v, T) == 1)
        return v;
    return boost::none;
  };

  bool a_ok, b_ok;
  for (auto l : info.leaves()) {
    auto next = [&](uint x) {return info.parent(x, l);};
    auto bl = info.branching(l);
    for (auto x : info.support(l))
      for (auto a = x, b = next(a), c = next(b); a != bl; a = b, b = c, c = next(c))
        if ((a_ok = a == x || is_branching(a)) || is_branching(c))
          if ((b_ok = is_branching(b)) || supported_other(b, l, x)) {
            add_edge(l, x, T);
            remove_edge(a_ok ? a : c, b, T);
            info.update();
            if (!b_ok) rule1action(b, *supported_other(b, l, x), T, info);
            return true;
          }
  }
  return false;
}

template<class Graph, class Tree, class LeafInfo>
bool rule19(Graph& G, Tree& T, LeafInfo& info) {
  // AKA rule16 v2
  auto is_branching = [&](uint x) {return out_degree(x, T) > 2;};
  std::tuple<uint, uint, uint> typedef Extra;
  vector<optional<Extra>> extra(num_vertices(G));

  for (auto l : info.leaves()) {
    auto next = [&](uint x) {return info.parent(x, l);};
    auto bl = info.branching(l);
    for (auto x : info.support(l)) {
      for (uint a = x, b = next(a), c = next(b); b != bl; a = b, b = c, c = next(c)) {
        if (!is_branching(b)) {
          if (x == a || is_branching(a)) {
            extra[b] = Extra(l, x, a);
          } else if (is_branching(c)) {
            extra[b] = Extra(l, x, c);
          }
        }
      }
    }
  }

  for (auto l2 : info.leaves()) {
    for (auto b : range(adjacent_vertices(l2, G))) {
      auto e = extra[b];
      if(e && get<0>(*e) != l2 && get<1>(*e) != l2) {
        uint l1, x, a;
        std::tie(l1, x, a) = *e;
        add_edge(l1, x, T);
        remove_edge(a, b, T);
        info.update();
        rule1action(b, l2, T, info);
        return true;
      }
    }
  }
  return false;
}

template<class Graph, class Tree>
std::function<vector<unsigned>(Graph&, Tree&)> make_improvement(std::string name) {
  vector<std::function<bool(Graph&, Tree&, leaf_info<Graph, Tree>&)>> typedef Rules;
  Rules rules = {
      rule0<Graph, Tree, leaf_info<Graph, Tree>>,
      rule1<Graph, Tree, leaf_info<Graph, Tree>>,
      rule2<Graph, Tree, leaf_info<Graph, Tree>>,
      rule3<Graph, Tree, leaf_info<Graph, Tree>>,
      rule4<Graph, Tree, leaf_info<Graph, Tree>>,
      rule5<Graph, Tree, leaf_info<Graph, Tree>>,
      rule6<Graph, Tree, leaf_info<Graph, Tree>>,
      rule7<Graph, Tree, leaf_info<Graph, Tree>>,
      rule8<Graph, Tree, leaf_info<Graph, Tree>>,
      rule9<Graph, Tree, leaf_info<Graph, Tree>>,
      rule10<Graph, Tree, leaf_info<Graph, Tree>>,
      rule11<Graph, Tree, leaf_info<Graph, Tree>>,
      rule12<Graph, Tree, leaf_info<Graph, Tree>>,
      rule13<Graph, Tree, leaf_info<Graph, Tree>>,
      rule14<Graph, Tree, leaf_info<Graph, Tree>>,
      rule15<Graph, Tree, leaf_info<Graph, Tree>>,
      rule16<Graph, Tree, leaf_info<Graph, Tree>>,
      rule17<Graph, Tree, leaf_info<Graph, Tree>>,
      rule18<Graph, Tree, leaf_info<Graph, Tree>>,
      rule19<Graph, Tree, leaf_info<Graph, Tree>>,
      rule20<Graph, Tree, leaf_info<Graph, Tree>>,
      rule21<Graph, Tree, leaf_info<Graph, Tree>>,
  };

  vector<unsigned> active;
  bool lazy = !found("eager", name);
  bool rand = found("rand", name);
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
    return [](Graph& G, Tree& T) {
      return vector<unsigned>(1, 0);
    };
  }
  else {
    std::stringstream ss(name);
    unsigned i;
    while (ss >> i) {
      active.push_back(i);
      if (ss.peek() == '+') ss.ignore();
    }
    if (active.empty()) throw std::invalid_argument("Unknown improvement method: " + name);
  }

  bool leafish = false;
  for (auto rule : active)
    if (9 > rule && rule < 15) leafish = true;

  return [=](Graph& G, Tree& T) mutable {
    leaf_info<Graph, Tree> info(G, T, leafish, lazy);
    bool applied = true;
    vector<unsigned> counter(active.size(), 0);
    auto order = active;
#ifndef NDEBUG
    unsigned steps = 0;
    //show("tree-" + std::to_string(steps) + ".dot", G, T);
#endif
    while (applied && !info.is_path()) {
      assert(++steps < num_vertices(G));
      assert(num_edges(T) == num_vertices(T) - 1);
      assert(is_connected(T));
      applied = false;
      if (rand) std::random_shuffle(order.begin(), order.end());
      for (unsigned i : order) {
        if (rules[i](G, T, info)) {
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

