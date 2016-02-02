// (C) 2014 Arek Olek

#pragma once

#include <iostream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include "graph.hpp"

struct graph_writer {
  bool labels;
  graph_writer(bool l) : labels(l) {}
  void operator()(std::ostream& out) const {
    if (labels) {
      out << "node [shape=circle color=burlywood style=filled]" << std::endl;
      out << "edge [style=dashed color=burlywood penwidth=4 weight=0.3]" << std::endl;
    } else {
      out << "node [shape=point color=burlywood style=filled label=\"\" width=0.05]" << std::endl;
      out << "edge [style=dashed color=burlywood]" << std::endl;
    }
  }
};

template<class Graph, class Tree>
class node_writer {
 public:
  node_writer(Graph const& g, Tree const & t)
      : T(t) {
    if (num_vertices(g) < 30 && is_planar(g)) {
      coords = straight_line_drawing(g);
    }
  }
  template<class V>
  void operator()(std::ostream& out, const V& v) const {
    out << "[";
    if (out_degree(v, T) == 1)
      out << "shape=square color=forestgreen";
    if (out_degree(v, T) > 2)
      out << "color=brown1";
    if (!coords.empty())
      out << " pos=\"" << coords[v].x << "," << coords[v].y << "!\"";
    out << "]";
  }
 private:
  Tree const& T;
  std::vector<coord_t> coords;
};

template<class Graph, class Tree>
node_writer<Graph, Tree> make_node_writer(Graph const & g, Tree const & t) {
  return node_writer<Graph, Tree>(g, t);
}

template<class Graph, class Tree>
class edge_writer {
public:
  edge_writer(Graph const& g, Tree const& t) : G(g), T(t) {}
  template <class E>
  void operator()(std::ostream& out, const E& e) const {
    if (edge(source(e, G), target(e, G), T).second)
      out << "[style=solid]";
  }
 private:
  Graph const& G;
  Tree const& T;
};

template<class Graph, class Tree>
edge_writer<Graph, Tree> make_edge_writer(Graph const & g, Tree const & t) {
  return edge_writer<Graph, Tree>(g, t);
}

template<class Graph, class Tree>
void show(std::string file, Graph const & g, Tree const & t) {
  std::ofstream f(file);
  write_graphviz(f, g, make_node_writer(g, t), make_edge_writer(g, t),
                 graph_writer(num_vertices(g) < 30));
  f.close();
}

template<class Graph>
void show(std::string file, Graph const & g) {
  std::ofstream f(file);
  write_graphviz(f, g, boost::default_writer(), boost::default_writer(),
                 graph_writer(false));
  f.close();
}
