// (C) 2014 Arek Olek

#pragma once

#include <algorithm>
#include <functional>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/function_input_iterator.hpp>
#include <boost/function_output_iterator.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Point_set_2.h>

#include "range.hpp"

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
void add_edges_uniform(Graph& G, double p, Generator generator) {
  std::bernoulli_distribution distribution(p);
  auto trial = std::bind(distribution, generator);
  unsigned n = num_vertices(G);
  for(unsigned i = 0; i < n; ++i)
    for(unsigned j = i + 1; j < n; ++j)
      if(trial())
        add_edge(i, j, G);
}

/*
def geometric(g, p, d=1.0):
    h = Graph(directed=False)
    h.add_vertex(g.num_vertices())
    def neighbors(s, vis):
        def visit(v):
            vis[v] = True
            for u in v.all_neighbours():
                if not vis[u] and np.linalg.norm(p[u].a-p[s].a) <= d:
                    if g.vertex_index[s] < g.vertex_index[u]:
                        h.add_edge(s, u)
                    visit(u)
        visit(s)
    for v in g.vertices():
        neighbors(v, g.new_vertex_property('bool'))
    return h
*/

template<class Graph, class Generator>
void add_random_geometric(Graph& G, double d, Generator generator) {
  CGAL::Exact_predicates_inexact_constructions_kernel           typedef Kernel;
  CGAL::Triangulation_vertex_base_with_info_2<unsigned, Kernel> typedef Vertex_struct;
  CGAL::Triangulation_data_structure_2<Vertex_struct>           typedef Triangulation_struct;
  CGAL::Point_set_2<Kernel, Triangulation_struct>               typedef Delaunay;
  CGAL::Circle_2<Kernel>                                        typedef Circle;
  Delaunay::Vertex_handle                                       typedef Vertex_handle;
  Delaunay::Point                                               typedef Point;

  std::uniform_real_distribution<> distribution(0, 1);
  auto r = std::bind(distribution, generator);

  //auto w = [&](Vertex_handle v){ std::cerr << v << std::endl; };

  int n = boost::num_vertices(G), counter = 0;
  std::function<std::pair<Point, int>()> f =
    [&]() {return std::make_pair(Point(r(), r()), counter++); };

  Delaunay D(boost::make_function_input_iterator(f, 0),
             boost::make_function_input_iterator(f, n));
  for(auto v : range(D.finite_vertices_begin(), D.finite_vertices_end())) {
    std::cerr << v.info() << " " << v << std::endl;
    Circle c(v.point(), d*d);
    //D.range_search(c, boost::make_function_output_iterator(w));
    std::vector<Vertex_handle> V;
    D.range_search(c, std::back_inserter(V));
    //for(auto x : V) std::cerr << x.info() << " ";
    //std::cerr << std::endl;
  }
}
