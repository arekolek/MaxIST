
#pragma once

#include <iostream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

using namespace std;

struct red_points_dashed {
  void operator()(ostream& out) const {
    out << "node [shape=point color=brown]" << endl;
    out << "edge [style=dashed color=burlywood1]" << endl;
  }
};

template <class Graph>
class node_writer {
public:
  node_writer(Graph const& t) : T(t) {}
  template <class V>
  void operator()(ostream& out, const V& v) const {
    if(degree(v, T) == 1)
      out << "[color=green]";
    if(degree(v, T) > 2)
      out << "[color=red]";
  }
private:
  Graph const& T;
};

template <class Graph>
node_writer<Graph> make_node_writer(Graph t) {
  return node_writer<Graph>(t);
}

template <class Graph>
class edge_writer {
public:
  edge_writer(Graph const& g, Graph const& t) : G(g), T(t) {}
  template <class E>
  void operator()(ostream& out, const E& e) const {
    auto vs = adjacent_vertices(source(e, G), T);
    for(auto vit = vs.first; vit != vs.second; ++vit)
      if(*vit == target(e, G)) {
        out << "[style=solid color=burlywood]";
        break;
      }
  }
private:
  Graph const& G;
  Graph const& T;
};

template <class Graph>
edge_writer<Graph> make_edge_writer(Graph g, Graph t) {
  return edge_writer<Graph>(g, t);
}

template<class Graph>
void show(string file, Graph g, Graph t) {
  ofstream f(file);
  write_graphviz(f, g,
    make_node_writer(t), make_edge_writer(g, t), red_points_dashed());
  f.close();
}

