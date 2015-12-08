// (C) 2014 Arek Olek

#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

#include <omp.h>

#include "bfs.hpp"
#include "dfs.hpp"
#include "fifodfs.hpp"
#include "fivethree.hpp"
#include "lost.hpp"
#include "rdfs.hpp"
#include "random.hpp"
#include "ilst.hpp"
#include "greedy.hpp"

#include "test_suite.hpp"

#include "debug.hpp"
#include "graph.hpp"
#include "options.hpp"
#include "range.hpp"
#include "timing.hpp"

boost::adjacency_matrix<boost::undirectedS> typedef amatrix;
boost::adjacency_list<boost::hash_setS, boost::vecS, boost::undirectedS> typedef ahash;
boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS> typedef alist;
boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> typedef aset;
boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> typedef avec;

template<class Graph, class Tree>
std::function<Tree(Graph&)> make_construction(std::string name) {
  if (name == "bfs")
    return bfs_tree<Graph, Tree> ;
  if (name == "dfs")
    return dfs_tree<Graph, Tree> ;
  if (name == "rdfs")
    return rdfs_tree<Graph, Tree> ;
  if (name == "rdfs50")
    return rdfs_best_tree<Graph, Tree> ;
  if (name == "fifo")
    return fifo_dfs_tree<Graph, Tree> ;
  if (name == "random")
    return random_tree<Graph, Tree> ;
  if (name == "wilson")
    return wilson_tree<Graph, Tree> ;
  if (name == "greedy")
    return greedy_tree<Graph, Tree> ;
  if (name == "ilst")
    return ilst<Graph, Tree> ;
  if (name == "5/3")
    return five_three_tree<Graph, Tree> ;
  throw std::invalid_argument("Unknown construction method: " + name);
}

template<class Graph>
std::function<unsigned(Graph&)> make_upper(std::string name) {
  if (name == "5/3") {
    return [](const Graph& G) {
      return 5*(num_vertices(G)-3)/6;
    };
  }
  return upper_bound<Graph> ;
}

template <class Graph, class Tree, class Suite, class Strings>
void run(Suite& suite, Strings const & constructions, Strings const & improvements, bool scratch) {
  #pragma omp parallel
  {
    auto id = omp_get_thread_num();
    double elapsed_c, elapsed_i;
    timing timer;
    #pragma omp for schedule(dynamic) nowait
    for(uint i = 0; i < suite.size(); ++i) {
      std::stringstream buffer;
      auto trial = suite.get(i);
      auto G = std::get<0>(trial);
      if(!is_connected(G)) continue;

      for(auto cname : constructions) {
        auto construct = make_construction<Graph, Tree>(cname);
        auto upper = make_upper<Graph>(cname);
        timer.start();
        auto T = construct(G);
        elapsed_c = timer.stop();

        assert(num_edges(T) == num_vertices(T)-1);

        const auto tree(T);

        for(auto iname : improvements) {
          if(scratch) T = tree;

          auto improve = make_improvement<Graph, Tree>(iname);
          timer.start();
          auto rules = improve(G, T);
          elapsed_i = timer.stop();

//  1    2       3      4       5      6         7      8      9             10           11        12     13     14     15
//  run  thread  model  degree  param  vertices  edges  upper  construction  improvement  internal  ctime  itime  steps  rules
          buffer
            << i << '\t'
            << id << '\t'
            << suite.type() << '\t'
            << std::get<1>(trial) << '\t'
            << std::get<2>(trial) << '\t'
            << num_vertices(G) << '\t'
            << num_edges(G) << '\t'
            << upper(G) << '\t'
            << cname << '\t'
            << iname << '\t'
            << num_internal(T) << '\t'
            << elapsed_c << '\t'
            << elapsed_i << '\t'
            << sum(rules) << '\t'
            << join(rules, "-") << '\t'
            << std::endl;
            ;
          //show("graph" + std::to_string(i) + ".dot", G, T);

          if(!scratch) elapsed_c += elapsed_i;
        }
      }
      #pragma omp critical
      std::cout << buffer.str() << std::flush;
    }
  }
}

template <class Graph, class Tree, class Sizes, class Params, class Strings>
void run(std::string t, unsigned z, Sizes sizes, Params params,
         Strings const & constructions, Strings const & improvements, bool scratch) {
  if (t.find(".xml") != std::string::npos) {
    //real_suite<Graph> suite(t, z);
    //run<Graph, Tree>(suite, constructions, improvements, scratch);
  } else if (t.find('.') != std::string::npos) {
    //file_suite<Graph> suite(t);
    //run<Graph, Tree>(suite, constructions, improvements, scratch);
  } else {
    for (auto n : sizes) {
      for (auto p : params) {
        test_suite<Graph> suite(t, z, n, p);
        run<Graph, Tree>(suite, constructions, improvements, scratch);
      }
    }
  }
}

int main(int argc, char** argv){
  using std::string;
  using std::vector;

  std::ios_base::sync_with_stdio(0);

  options opt(argc, argv);
  auto z = opt.get<int>("-z", 1);
  auto sizes = opt.getList<int>("-n", {100});
  auto tests = opt.getList<string>("-t", {"gnp+mst"});
  auto degrees = opt.getList<float>("-d", { 3. });
  auto constructions = opt.getList<string>("-c", {"bfs", "dfs", "rdfs", "fifo", "rdfs50", "ilst", "random"});
  auto improvements = opt.getList<string>("-i", {"none", "prieto", "lost-light", "lost", "lost-ex"});
  auto scratch = opt.has("--scratch");

  for(auto t : tests)
    run<alist, alist>(t, z, sizes, degrees, constructions, improvements, scratch);

  return 0;
}
