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
#include "random.hpp"

#include "test_suite.hpp"

#include "debug.hpp"
#include "graph.hpp"
#include "options.hpp"
#include "range.hpp"
#include "timing.hpp"

boost::property<boost::edge_color_t, boost::default_color_type> typedef color;
boost::adjacency_list<boost::hash_setS, boost::vecS, boost::undirectedS,
    boost::no_property, color> typedef graph;

template<class Graph>
std::function<Graph(Graph&)> make_construction(std::string name) {
  if (name == "bfs")
    return bfs_tree<Graph> ;
  if (name == "dfs")
    return dfs_tree<Graph> ;
  if (name == "rdfs")
    return rdfs_tree<Graph> ;
  if (name == "rdfs500")
    return rdfs_best_tree<Graph> ;
  if (name == "fifo")
    return fifo_dfs_tree<Graph> ;
  if (name == "random")
      return random_tree<Graph> ;
  if (name == "5/3")
    return five_three_tree<Graph> ;
  throw std::invalid_argument("Unknown construction method: " + name);
}

template <class Graph, class Suite>
void run(Suite& suite, std::string cname, std::string iname) {
  auto construct = make_construction<Graph>(cname);
  auto improve = make_improvement<Graph,Graph>(iname);
  #pragma omp parallel
  {
    auto id = omp_get_thread_num();
    std::stringstream buffer;
    unsigned steps;
    double elapsed_c, elapsed_i;
    timing timer;
    #pragma omp for schedule(dynamic) nowait
    for(uint i = 0; i < suite.size(); ++i) {
        auto G = suite.get(i);
        if(!is_connected(G)) {
          continue;
        }
        timer.start();
        auto T = construct(G);
        elapsed_c = timer.stop();
        assert(num_edges(T) == num_vertices(T)-1);
        timer.start();
        steps = improve(G, T);
        elapsed_i = timer.stop();
        // run type parameter vertices edges upper construction improvement internal time steps
        buffer
          << i << '\t'
          << id << '\t'
          << suite.type() << '\t'
          << suite.parameter() << '\t'
          << num_vertices(G) << '\t'
          << num_edges(G) << '\t'
          << upper_bound(G) << '\t'
          << cname << '\t'
          << iname << '\t'
          << num_internal(T) << '\t'
          << elapsed_c << '\t'
          << elapsed_i << '\t'
          << steps << '\t'
          << std::endl;
          ;
        //show("graph" + to_string(i) + ".dot", G, T);
    }
    #pragma omp critical
    std::cout << buffer.str();
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

  //std::cout << "run,type,parameter,vertices,edges,upper,construction,improvement,internal,time,steps" << std::endl;
  for(auto t : tests)
    for(auto c : constructions)
      for(auto i : improvements)
        run<graph>(t, z, sizes, parameters, c, i);

  return 0;
}
