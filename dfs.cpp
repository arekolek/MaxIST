// (C) 2014 Arek Olek

#include <iostream>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/undirected_dfs.hpp>

#include "debug.hpp"

using std::cin;
using std::cout;
using std::endl;
using std::vector;

using boost::edge_color;
using boost::on_tree_edge;

boost::property<boost::edge_color_t,
  boost::default_color_type>            typedef color;
boost::adjacency_list<
  boost::hash_setS, boost::vecS, boost::undirectedS,
  boost::no_property, color>            typedef graph;

template <class ParentDecorator>
struct print_parent {
  print_parent(const ParentDecorator& p_) : p(p_) { }
  template <class Vertex>
  void operator()(const Vertex& v) const {
    std::cout << "parent[" << v << "] = " <<  p[v]  << std::endl;
  }
  ParentDecorator p;
};


template <class NewGraph, class Tag>
struct graph_copier
  : public boost::base_visitor<graph_copier<NewGraph, Tag> >
{
  typedef Tag event_filter;

  graph_copier(NewGraph& graph) : new_g(graph) { }

  template <class Edge, class Graph>
  void operator()(Edge e, Graph& g) {
    boost::add_edge(boost::source(e, g), boost::target(e, g), new_g);
  }
private:
  NewGraph& new_g;
};

template <class NewGraph, class Tag>
inline graph_copier<NewGraph, Tag>
copy_graph(NewGraph& g, Tag) {
  return graph_copier<NewGraph, Tag>(g);
}

int main(int argc, char** argv){
  std::ios_base::sync_with_stdio(0);

  int Z, sum = 0;
  cin >> Z;

  for(int z = 0; z < Z; ++z) {
    int n, m;
    cin >> n >> m;

    graph G(n);
    for(int i = 0; i < m; ++i) {
      unsigned a, b;
      cin >> a >> b;
      add_edge(a, b, G);
    }

    vector<unsigned> p(num_vertices(G));
    for(unsigned v = 0; v < p.size(); ++v) p[v] = v;

    graph T;

    undirected_dfs(G,
      visitor(make_dfs_visitor(copy_graph(T, on_tree_edge())))
        .edge_color_map(get(edge_color, G)));

    int internal = 0;
    int upper = n - 2;
    auto vs = vertices(T);
    for(auto vit = vs.first; vit != vs.second; ++vit)
      internal += boost::degree(*vit, T) > 1;

    cout << n-internal << " " << internal << " " << upper
      << " " << n << " " << (upper-internal)/(double)upper << endl;

    sum += n-internal;

    show("graph.dot", G, T);

  }

  cout << sum/(double)Z << endl;

  return 0;
}

