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
boost::adjacency_list<boost::slistS, boost::vecS, boost::undirectedS> typedef aslist;
boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> typedef aset;
boost::adjacency_list<boost::multisetS, boost::vecS, boost::undirectedS> typedef amultiset;
boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> typedef avec;

template<class Graph, class Tree>
std::function<Tree(Graph&)> make_construction(std::string name,
                                              unsigned size,
                                              unsigned index,
                                              std::string seed) {
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
  if (name == "5/3") {
    return [=](const Graph& G) {
      std::default_random_engine generator(index);
      Tree g(num_vertices(G));
      copy_edges_shuffled(five_three_tree<Graph, Tree>(G), g, generator);
      return g;
    };
  }
  if (found(".xml", name)) {
    return [=](const Graph& G) {
      real_suite<Tree> suite(name, size, seed);
      return std::get<0>(suite.get(index));
    };
  }
  throw std::invalid_argument("Unknown construction method: " + name);
}

template<class Graph>
std::function<unsigned(Graph&)> make_upper(std::string name) {
  if (name == "53") {
    return [](const Graph& G) {
      return 5*(num_vertices(G)-3)/6;
    };
  }
  return upper_bound<Graph> ;
}

std::string type(amatrix const & G) { return "matrix"; }
std::string type(ahash const & G) { return "unordered_set"; }
std::string type(aset const & G) { return "set"; }
std::string type(amultiset const & G) { return "multiset"; }
std::string type(alist const & G) { return "list"; }
std::string type(aslist const & G) { return "slist"; }
std::string type(avec const & G) { return "vector"; }

template <class Graph, class Tree, class Suite, class Strings>
void run(Suite& suite,
         Strings const & constructions,
         Strings const & improvements,
         bool scratch,
         std::string seed) {
  timing total(CLOCK_REALTIME);
  total.start();
  auto graph_type = type(Graph(0));
  auto tree_type = type(Tree(0));
  #pragma omp parallel
  {
    auto id = omp_get_thread_num();
    double elapsed_c, elapsed_i;
    std::stringstream buffer;
    timing timer;
    #pragma omp for schedule(dynamic) nowait
    for(uint i = 0; i < suite.size(); ++i) {
      auto trial = suite.get(i);
      auto G = std::get<0>(trial);

      for (auto cname : constructions) {
        auto construct = make_construction<Graph, Tree>(cname, suite.size(), suite.get_seed(i), seed);
        auto upper = make_upper<Graph>(suite.type());
        timer.start();
        auto T = construct(G);
        elapsed_c = timer.stop();

        assert(num_edges(T) == num_vertices(T) - 1);
        assert(is_connected(T));
        assert(is_subgraph(T, G));

        const auto tree(T);

        for (auto iname : improvements) {
          if (scratch) T = tree;
          auto num_internal_before = num_internal(T);
          auto improve = make_improvement<Graph, Tree>(iname);

          timer.start();
          auto rules = improve(G, T);
          elapsed_i = timer.stop();

          buffer
#ifdef GRAPH_REPR
            << graph_type << '\t'
#endif
#ifdef TREE_REPR
            << tree_type << '\t'
#endif
            << std::get<1>(trial) << '\t'
            << id << '\t'
            << suite.type() << '\t'
            << std::get<2>(trial) << '\t'
            << std::get<3>(trial) << '\t'
            << num_vertices(G) << '\t'
            << num_edges(G) << '\t'
            << upper(G) << '\t'
            << cname << '\t'
            << iname << '\t'
            << num_internal_before << '\t'
            << num_internal(T) << '\t'
            << elapsed_c << '\t'
            << elapsed_i << '\t'
            << sum(rules) << '\t'
            << join(rules, "-")
            << std::endl;
            ;
#ifndef NDEBUG
          if (suite.size() < 10)
            show("graph" + std::to_string(i) + graph_type + ".dot", G, T);
          //if(rules[rules.size()-2] > 10) show("graph" + std::to_string(i) + ".dot", G, T);
#endif
          if (!scratch) elapsed_c += elapsed_i;
        }
      }

      if (id == 0 || i + 1 == suite.size()) {
        std::cerr << '\t' << 100 * (i + 1) / suite.size() << "% of "
        << suite.size();
        std::cerr << '\t' << readable(total.stop()) << " elapsed    \r";
      }
    }
    #pragma omp critical
    std::cout << buffer.str() << std::flush;
  }
  std::cerr << std::endl;
}

template <class Graph, class Tree, class Sizes, class Params, class Strings>
void run(std::string t,
         unsigned z,
         Sizes sizes,
         Params params,
         Strings const & constructions,
         Strings const & improvements,
         bool scratch,
         std::string seed) {
  if (found(".xml", t)) {
    real_suite<Graph> suite(t, z, seed);
    run<Graph, Tree>(suite, constructions, improvements, scratch, seed);
  } else if (found(".", t)) {
    file_suite<Graph> suite(t, z, seed);
    run<Graph, Tree>(suite, constructions, improvements, scratch, seed);
  } else {
    test_suite<Graph> suite(t, z, sizes, params, seed);
    run<Graph, Tree>(suite, constructions, improvements, scratch, seed);
  }
}

template <class Graph, class Sizes, class Degrees, class Strings>
void run(std::string tree,
         std::string model,
         unsigned sample,
         Sizes sizes,
         Degrees degrees,
         Strings const & constructions,
         Strings const & improvements,
         bool scratch,
         std::string seed) {
#ifdef TREE_REPR
  if (tree == "matrix")
    run<Graph, amatrix>  (model, sample, sizes, degrees, constructions, improvements, scratch, seed);
  if (tree == "unordered_set")
    run<Graph, ahash>    (model, sample, sizes, degrees, constructions, improvements, scratch, seed);
  if (tree == "multiset")
    run<Graph, amultiset>(model, sample, sizes, degrees, constructions, improvements, scratch, seed);
  if (tree == "list")
    run<Graph, alist>    (model, sample, sizes, degrees, constructions, improvements, scratch, seed);
  if (tree == "set")
    run<Graph, aset>     (model, sample, sizes, degrees, constructions, improvements, scratch, seed);
  if (tree == "vector")
    run<Graph, avec>     (model, sample, sizes, degrees, constructions, improvements, scratch, seed);
#endif
  if (tree == "slist")
    run<Graph, aslist>   (model, sample, sizes, degrees, constructions, improvements, scratch, seed);
}

template <class Sizes, class Degrees, class Strings>
void run(std::string graph,
         std::string tree,
         std::string model,
         unsigned sample,
         Sizes sizes,
         Degrees degrees,
         Strings const & constructions,
         Strings const & improvements,
         bool scratch,
         std::string seed) {
#ifdef GRAPH_REPR
  if (graph == "matrix")
    run<amatrix>  (tree, model, sample, sizes, degrees, constructions, improvements, scratch, seed);
  if (graph == "unordered_set")
    run<ahash>    (tree, model, sample, sizes, degrees, constructions, improvements, scratch, seed);
  if (graph == "multiset")
    run<amultiset>(tree, model, sample, sizes, degrees, constructions, improvements, scratch, seed);
  if (graph == "slist")
    run<aslist>   (tree, model, sample, sizes, degrees, constructions, improvements, scratch, seed);
  if (graph == "vector")
    run<avec>     (tree, model, sample, sizes, degrees, constructions, improvements, scratch, seed);
  if (graph == "list")
    run<alist>    (tree, model, sample, sizes, degrees, constructions, improvements, scratch, seed);
#endif
  if (graph == "set")
    run<aset>     (tree, model, sample, sizes, degrees, constructions, improvements, scratch, seed);
}

std::string time_seed() {
  return std::to_string(std::chrono::system_clock::now().time_since_epoch().count());
}

int main(int argc, char** argv){
  using std::string;
  using std::vector;

  std::ios_base::sync_with_stdio(0);

  options opt(argc, argv);
  auto z = opt.get<int>("-z", 100);
  auto graphs = opt.getList<string>("-g", {"set"});
  auto trees = opt.getList<string>("-t", {"slist"});
  auto sizes = opt.getList<unsigned>("-n", {1000});
  auto tests = opt.getList<string>("-m", {"gnp+mst"});
  auto degrees = opt.getList<double>("-d", { 3. });
  auto constructions = opt.getList<string>("-c", {"wilson"});
  auto improvements = opt.getList<string>("-i", {"none"});
  auto scratch = opt.has("--scratch");
  auto seed = opt.get<string>("-s", time_seed());

  std::cout << "# Seed: " << seed << std::endl;

  std::cout
#ifdef GRAPH_REPR
    << "graph\t"
#endif
#ifdef TREE_REPR
    << "tree\t"
#endif
    << "run\t"
    << "thread\t"
    << "model\t"
    << "degree\t"
    << "param\t"
    << "vertices\t"
    << "edges\t"
    << "upper\t"
    << "construction\t"
    << "improvement\t"
    << "before\t"
    << "internal\t"
    << "ctime\t"
    << "itime\t"
    << "steps\t"
    << "rules"
    << std::endl;

  for(auto t : tests)
    for(auto graph : graphs)
      for(auto tree : trees)
        run(graph, tree, t, z, sizes, degrees, constructions, improvements,
            scratch, seed);

  return 0;
}
