//=======================================================================
// Copyright 2007 Aaron Windsor
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
#include <iostream>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <boost/graph/planar_canonical_ordering.hpp>
#include <boost/graph/is_straight_line_drawing.hpp>
#include <boost/graph/chrobak_payne_drawing.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/make_connected.hpp>
#include <boost/graph/make_biconnected_planar.hpp>
#include <boost/graph/make_maximal_planar.hpp>

#include "graph.hpp"
#include "range.hpp"

//a class to hold the coordinates of the straight line embedding
struct coord_t {
  std::size_t x;
  std::size_t y;
};

template<class G>
std::vector<coord_t> straight_line_drawing(G const & gIn) {
  using namespace boost;

  typedef adjacency_list<vecS, vecS, undirectedS, property<vertex_index_t, int>,
      property<edge_index_t, int> > Graph;

  Graph g;
  copy_edges(gIn, g);

  //Define the storage type for the planar embedding
  typedef std::vector<std::vector<graph_traits<Graph>::edge_descriptor>> embedding_storage_t;
  typedef iterator_property_map<embedding_storage_t::iterator,
      property_map<Graph, vertex_index_t>::type> embedding_t;

  // Create the planar embedding
  embedding_storage_t embedding_storage(num_vertices(g));
  embedding_t embedding(embedding_storage.begin(), get(vertex_index, g));

  boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                               boyer_myrvold_params::embedding = embedding);

  make_connected(g);
  make_biconnected_planar(g, embedding);
  make_maximal_planar(g, embedding);

  // Find a canonical ordering
  std::vector<typename graph_traits<Graph>::vertex_descriptor> ordering;
  planar_canonical_ordering(g, embedding, std::back_inserter(ordering));

  //Set up a property map to hold the mapping from vertices to coord_t's
  typedef std::vector<coord_t> drawing_storage_t;
  typedef boost::iterator_property_map<drawing_storage_t::iterator,
      property_map<Graph, vertex_index_t>::type> drawing_t;

  drawing_storage_t drawing_storage(num_vertices(g));
  drawing_t drawing(drawing_storage.begin(), get(vertex_index, g));

  // Compute the straight line drawing
  chrobak_payne_straight_line_drawing(g, embedding, ordering.begin(),
                                      ordering.end(), drawing);

  return drawing_storage;
}

int main(int argc, char** argv) {
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph;

  // you can use the functions make_connected, make_biconnected_planar,
  // and make_maximal planar in sequence to add a set of edges
  // to any undirected planar graph to make it maximal planar.

  graph g(7);
  add_edge(0, 1, g);
  add_edge(1, 2, g);
  add_edge(2, 3, g);
  add_edge(3, 0, g);
  //add_edge(3, 4, g);
  add_edge(4, 5, g);
  add_edge(5, 6, g);
  //add_edge(6, 3, g);
  add_edge(0, 4, g);
  add_edge(1, 3, g);
  add_edge(3, 5, g);
  add_edge(2, 6, g);
  add_edge(1, 4, g);
  add_edge(1, 5, g);
  add_edge(1, 6, g);

  auto coords = straight_line_drawing(g);

  std::cout << "The straight line drawing is: " << std::endl;
  for (auto v : ::range(vertices(g))) {
    coord_t coord(coords[v]);
    std::cout << v << " -> (" << coord.x << ", " << coord.y << ")" << std::endl;
  }

  return 0;
}
