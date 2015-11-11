// (C) 2014 Arek Olek

#pragma once

#include <algorithm>
#include <functional>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/function_input_iterator.hpp>
#include <boost/function_output_iterator.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Point_set_2.h>
#include <CGAL/squared_distance_2.h>

#include "range.hpp"

template<class Graph>
unsigned num_internal(Graph const & G) {
  unsigned internal = 0;
  for (auto v : range(vertices(G)))
    internal += out_degree(v, G) > 1;
  return internal;
}

template<class Graph>
unsigned upper_bound(Graph const & G) {
  return std::min(num_internal(G), (unsigned)num_vertices(G) - 2);
}

template<class Graph>
bool is_connected(Graph const & G) {
  std::vector<int> component(num_vertices(G));
  return connected_components(G, &component[0]) == 1;
}

template <class Input, class Output>
void copy_edges(const Input& in, Output& out) {
  for(auto e : range(edges(in)))
    add_edge(source(e, in), target(e, in), out);
}

template<class Graph, class Generator>
void add_spider(Graph& G, unsigned legs, Generator generator) {
  unsigned n = num_vertices(G);
  unsigned cutoff = n / legs;
  auto range_iterator = boost::make_counting_iterator<int>;
  std::vector<int> path(range_iterator(0), range_iterator(n));
  std::shuffle(path.begin(), path.end(), generator);
  for(unsigned i = 0; i < n-1; ++i)
    add_edge(path[i % cutoff == 0 ? 0 : i], path[i+1], G);
}

template<class Graph, class Generator>
void add_edges_uniform(Graph& G, double p, Generator generator, bool mst) {
  unsigned n = num_vertices(G);
  double connectedness = 20.0 / n;
  if(mst && p < connectedness) {
    typedef boost::property<boost::edge_weight_t, double> Weight;
    typedef boost::adjacency_list<boost::hash_setS, boost::vecS, boost::undirectedS, boost::no_property, Weight> WeightedGraph;
    typedef boost::graph_traits<WeightedGraph>::edge_descriptor WeightedEdge;
    std::uniform_real_distribution<> distribution(0, 1);
    auto trial = std::bind(distribution, generator);
    WeightedGraph g;
    for(unsigned i = 0; i < n; ++i) {
      for(unsigned j = i + 1; j < n; ++j) {
        auto w = trial();
        if(w < connectedness) add_edge(i, j, w, g);
        if(w < p) add_edge(i, j, G);
      }
    }
    std::vector<WeightedEdge> t;
    kruskal_minimum_spanning_tree(g, std::back_inserter(t));
    for(auto e : t) add_edge(source(e, g), target(e, g), G);
  } else {
    std::bernoulli_distribution distribution(p);
    auto trial = std::bind(distribution, generator);
    for(unsigned i = 0; i < n; ++i) {
      for(unsigned j = i + 1; j < n; ++j) {
        if(trial()) add_edge(i, j, G);
      }
    }
  }
}

class Geometric {
  CGAL::Exact_predicates_inexact_constructions_kernel           typedef Kernel;
  CGAL::Triangulation_vertex_base_with_info_2<unsigned, Kernel> typedef Vertex_struct;
  CGAL::Triangulation_data_structure_2<Vertex_struct>           typedef Triangulation_struct;
  CGAL::Point_set_2<Kernel, Triangulation_struct>               typedef Delaunay;
  Delaunay::Point                                               typedef Point;

  Delaunay D;

public:
  template<class Generator>
  Geometric(int n, Generator generator) {
    // Choose random points in square [0, 1) x [0, 1)
    std::uniform_real_distribution<> distribution(0, 1);
    auto r = std::bind(distribution, generator);
    // Generate index for every point
    int counter = 0;
    std::function<std::pair<Point, int>()> f =
      [&]() {return std::make_pair(Point(r(), r()), counter++); };
    // Initialize triangulation with n points
    D = Delaunay(boost::make_function_input_iterator(f, 0),
                 boost::make_function_input_iterator(f, n));
  }

  template<class Graph>
  void add_random_geometric(Graph& G, double d) {
    CGAL::Circle_2<Kernel>  typedef Circle;
    Delaunay::Vertex_handle typedef Vertex_handle;
    // Query each point for other points lying in the circle
    for(auto v : range(D.finite_vertices_begin(), D.finite_vertices_end()))
      D.range_search(Circle(v.point(), d*d), boost::make_function_output_iterator(
          [&](Vertex_handle u){
            if(v.info() < u->info())
              add_edge(v.info(), u->info(), G);
      }));
  }

  template<class Graph>
  void add_mst(Graph& G) {
    typedef boost::property<boost::edge_weight_t, Kernel::FT> Weight;
    typedef boost::adjacency_list<boost::hash_setS, boost::vecS, boost::undirectedS, boost::no_property, Weight> WeightedGraph;
    typedef boost::graph_traits<WeightedGraph>::edge_descriptor WeightedEdge;

    Delaunay::Finite_vertices_iterator vit;
    Delaunay::Vertex_circulator vc, done;

    WeightedGraph g;

    for(vit = D.finite_vertices_begin(); vit != D.finite_vertices_end(); ++vit) {
      unsigned s = vit->info();
      done = vc = vit->incident_vertices();
      if(vc != 0) do {
        if(D.is_infinite(vc)) continue;
        unsigned d = vc->info();
        add_edge(s, d, CGAL::squared_distance(vit->point(), vc->point()), g);
      } while(++vc != done);
    }

    std::vector<WeightedEdge> mst;
    kruskal_minimum_spanning_tree(g, std::back_inserter(mst));
    for(auto e : mst) add_edge(source(e, g), target(e, g), G);
  }
};
