// (C) 2014 Arek Olek

#include <unordered_map>

#include <boost/graph/adjacency_matrix.hpp>

namespace detail {
  boost::adjacency_matrix<boost::undirectedS> typedef graph;
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
        std::cerr << "leaf " << l << std::endl;
        while(degree(x, T) == 2) {
          std::cerr << x << std::endl;
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
  std::cerr << "dfs" << std::endl;
  auto T = dfs_tree(G);
  std::cerr << "info" << std::endl;
  leaf_info<Graph> info(T);
  while(!info.is_path() && rule2(G, T, info));
  return T;
}

template <class Graph, class Tree, class LeafInfo>
bool rule2(Graph& G, Tree& T, LeafInfo& info) {
  std::cerr << "rule" << std::endl;
  for(auto x : info.leaves())
    for(auto y : info.leaves())
      if(edge(x, y, G).second) {
        // add x - y
        add_edge(x, y, T);
        // remove b(x) - b(x)->x
        auto bxx = info.branching_neighbor(x);
        remove_edge(info.branching(x), bxx, T);
        // update leaf list
        //info.remove_leaf(x);
        //info.remove_leaf(y);
        //if(degree(bxx, T) == 1)
          //info.add_leaf(bxx);
        std::cerr << "update" << std::endl;
        info.update();
        return true;
      }
  return false;
}
