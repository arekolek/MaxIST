// (C) 2014 Arek Olek

#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include <omp.h>

#include "bfs.hpp"
#include "dfs.hpp"
#include "fifodfs.hpp"
#include "fivethree.hpp"
#include "lost.hpp"
#include "rdfs.hpp"

#include "test_suite.hpp"

#include "debug.hpp"
#include "options.hpp"
#include "range.hpp"
#include "timing.hpp"

boost::property<boost::edge_color_t, boost::default_color_type> typedef color;
boost::adjacency_list<boost::hash_setS, boost::vecS, boost::undirectedS,
    boost::no_property, color> typedef graph;

template<class Graph>
unsigned num_internal(Graph const & G) {
  unsigned internal = 0;
  for (auto v : range(vertices(G)))
    internal += degree(v, G) > 1;
  return internal;
}

template<class Graph>
unsigned upper_bound(Graph const & G) {
  return std::min(num_internal(G), (unsigned)num_vertices(G) - 2);
}

template<class Graph>
std::function<Graph(Graph&)> make_construction(std::string name) {
  if (name == "bfs")
    return bfs_tree<Graph> ;
  if (name == "dfs")
    return dfs_tree<Graph> ;
  if (name == "rdfs")
    return rdfs_tree<Graph> ;
  if (name == "fifo")
    return fifo_dfs_tree<Graph> ;
  if (name == "rdfs-sort")
    return rdfs_sort_tree<Graph> ;
  if (name == "rdfs-rand")
    return rdfs_rand_tree<Graph> ;
  if (name == "5/3")
    return five_three_tree<Graph> ;
  throw std::invalid_argument("Unknown construction method: " + name);
}

template <class Graph, class Suite>
void run(Suite& suite, std::string cname, std::string iname) {
  auto construct = make_construction<Graph>(cname);
  auto improve = make_improvement<Graph,Graph>(iname);
  timing timer;
  unsigned steps;
  double elapsed;
  std::stringstream buf;
  #pragma omp parallel private(timer, steps, elapsed, buf)
  {
    #pragma omp for nowait
    for(uint i = 0; i < suite.size(); ++i) {
        auto G = suite.get(i);
        timer.start();
        auto T = construct(G);
        if(num_vertices(T) != num_vertices(G))
          // G is unconnected
          continue;
        // assert T is a spanning tree
        assert(num_edges(T) == num_vertices(T)-1);
        steps = improve(G, T);
        elapsed = timer.stop();
        // run type parameter vertices edges upper construction improvement internal time steps
        buf
          << i << ','
          << suite.type() << ','
          << suite.parameter() << ','
          << num_vertices(G) << ','
          << num_edges(G) << ','
          << upper_bound(G) << ','
          << cname << ','
          << iname << ','
          << num_internal(T) << ','
          << elapsed << ','
          << steps << std::endl;
          ;
        //show("graph" + to_string(i) + ".dot", G, T);
    }
    #pragma omp critical
    std::cout << buf.rdbuf();
  }
}

template <class Graph, class Sizes, class Params>
void run(std::string t, unsigned z, Sizes sizes, Params params, std::string cname, std::string iname) {
  if(t.find('.') != std::string::npos) {
    file_suite<Graph> suite(t);
    run<Graph>(suite, cname, iname);
  }
  else {
    for(auto n : sizes) {
      for(auto p : params) {
        test_suite<Graph> suite(t, z, n, p);
        run<Graph>(suite, cname, iname);
      }
    }
  }
}

int main(int argc, char** argv){
  using std::string;
  using std::vector;

  std::ios_base::sync_with_stdio(0);

  options opt(argc, argv);
  auto z = opt.get<int>("-z", 100);
  auto sizes = opt.getList<int>("-n", {100});
  auto tests = opt.getList<string>("-t", {"gnp"/*, "rgg"*/});
  auto parameters = opt.getList<float>("-p", {
      0.0001f, 0.0005f, 0.001f, 0.003f, 0.005f, 0.008f, 0.01f, 0.03f, 0.05f, 0.08f, 0.1f, 0.2f, 0.25f, 0.3f, 0.35f,
      //0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99
  });
  auto constructions = opt.getList<string>("-c", {"bfs", "dfs", "rdfs", "fifo", "rdfs-sort", "rdfs-rand"});
  auto improvements = opt.getList<string>("-i", {"none"/*, "prieto", "lost-light", "lost"*/});

  //~ vector<float> ps {0.0002, 0.0105, 0.021, 0.0312, 0.0415, 0.0518, 0.062, 0.0725, 0.0827};

  std::cout << "run,type,parameter,vertices,edges,upper,construction,improvement,internal,time,steps" << std::endl;
  for(auto t : tests)
    for(auto c : constructions)
      for(auto i : improvements)
        run<graph>(t, z, sizes, parameters, c, i);

  return 0;
}
