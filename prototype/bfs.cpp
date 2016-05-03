//=======================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================

#include <boost/config.hpp>

#include <algorithm>
#include <vector>
#include <utility>
#include <iostream>

#include <boost/graph/visitors.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graph_utility.hpp>

template<class NewGraph, class Tag>
struct graph_copier
: public boost::base_visitor<graph_copier<NewGraph, Tag> >
{
  typedef Tag event_filter;

  graph_copier(NewGraph& graph)
      : new_g(graph) {
  }

  template<class Edge, class Graph>
  void operator()(Edge e, Graph& g) {
    boost::add_edge(boost::source(e, g), boost::target(e, g), new_g);
  }
 private:
  NewGraph& new_g;
};

template<class NewGraph, class Tag>
inline graph_copier<NewGraph, Tag>
copy_graph(NewGraph& g, Tag) {
  return graph_copier<NewGraph, Tag>(g);
}

int main(int , char* [])
{
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

  Graph G(5);
  boost::add_edge(0, 2, G);
  boost::add_edge(1, 1, G);
  boost::add_edge(1, 3, G);
  boost::add_edge(1, 4, G);
  boost::add_edge(2, 1, G);
  boost::add_edge(2, 3, G);
  boost::add_edge(2, 4, G);
  boost::add_edge(3, 1, G);
  boost::add_edge(3, 4, G);
  boost::add_edge(4, 0, G);
  boost::add_edge(4, 1, G);

  Graph G_copy(5);


  // The source vertex
  boost::depth_first_search(G,
                            boost::visitor(boost::make_dfs_visitor(copy_graph(G_copy, boost::on_tree_edge()))));

  boost::print_graph(G);
  boost::print_graph(G_copy);

  return 0;
}
