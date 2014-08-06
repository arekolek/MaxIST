// (C) 2014 Arek Olek

#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/undirected_dfs.hpp>

using boost::edge_color;
using boost::on_tree_edge;

template <class NewGraph, class Tag>
struct graph_copier
  : public boost::base_visitor<graph_copier<NewGraph, Tag> >
{
  typedef Tag event_filter;

  graph_copier(NewGraph& graph) : new_g(graph) { }

  template <class Edge, class Graph>
  void operator()(Edge e, Graph& g) {
    add_edge(source(e, g), target(e, g), new_g);
  }
private:
  NewGraph& new_g;
};

template <class NewGraph, class Tag>
inline graph_copier<NewGraph, Tag>
copy_visitor(NewGraph& g, Tag) {
  return graph_copier<NewGraph, Tag>(g);
}

template <class Graph>
Graph dfs_tree(Graph& G) {
  Graph T;

  undirected_dfs(G,
    visitor(make_dfs_visitor(copy_visitor(T, on_tree_edge())))
      .edge_color_map(get(edge_color, G)));

  return T;
}
