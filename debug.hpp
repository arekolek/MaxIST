// (C) 2014 Arek Olek

#pragma once

#include <iostream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

struct red_points_dashed {
  void operator()(std::ostream& out) const {
    out << "node [shape=point color=brown]" << std::endl;
    out << "edge [style=dashed color=burlywood1]" << std::endl;
  }
};

template <class Graph>
class node_writer {
public:
  node_writer(Graph const& t) : T(t) {}
  template <class V>
  void operator()(std::ostream& out, const V& v) const {
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

template<class Graph, class Tree>
class edge_writer {
public:
  edge_writer(Graph const& g, Tree const& t) : G(g), T(t) {}
  template <class E>
  void operator()(std::ostream& out, const E& e) const {
    auto vs = adjacent_vertices(source(e, G), T);
    for(auto vit = vs.first; vit != vs.second; ++vit)
      if(*vit == target(e, G)) {
        out << "[style=solid color=burlywood]";
        break;
      }
  }
private:
  Graph const& G;
  Tree const& T;
};

template<class Graph, class Tree>
edge_writer<Graph, Tree> make_edge_writer(Graph g, Tree t) {
  return edge_writer<Graph, Tree>(g, t);
}

template<class Graph, class Tree>
void show(std::string file, Graph g, Tree t) {
  std::ofstream f(file);
  write_graphviz(f, g,
    make_node_writer(t), make_edge_writer(g, t), red_points_dashed());
  f.close();
}

