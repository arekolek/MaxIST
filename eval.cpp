// (C) 2014 Arek Olek

#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "bfs.hpp"
#include "dfs.hpp"
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
  if(name == "bfs")
    return bfs_tree<Graph>;
  if(name == "dfs")
      return dfs_tree<Graph>;
  if(name == "rdfs")
      return rdfs_tree<Graph>;
  throw std::invalid_argument("Unknown construction method: " + name);
}

template<class Graph, class Tree>
std::function<int(Graph&,Tree&)> make_improvement(std::string name) {
  std::vector<std::function<bool(Graph&,Tree&,leaf_info<Tree>&)>> typedef Rules;
  Rules rules = {
          rule1<Graph,Tree,leaf_info<Tree>>,
          rule2<Graph,Tree,leaf_info<Tree>>,
          rule3<Graph,Tree,leaf_info<Tree>>,
          rule4<Graph,Tree,leaf_info<Tree>>,
          rule5<Graph,Tree,leaf_info<Tree>>,
          rule6<Graph,Tree,leaf_info<Tree>>,
          rule7<Graph,Tree,leaf_info<Tree>>,
          rule8<Graph,Tree,leaf_info<Tree>>,
          rule9<Graph,Tree,leaf_info<Tree>>,
          rule10<Graph,Tree,leaf_info<Tree>>,
          rule11<Graph,Tree,leaf_info<Tree>>,
          rule12<Graph,Tree,leaf_info<Tree>>,
          rule13<Graph,Tree,leaf_info<Tree>>,
        };
  if (name == "prieto")
    rules = Rules(rules.begin() + 1, rules.begin() + 2);
  else if (name == "lost-light")
    rules = Rules(rules.begin() + 1, rules.begin() + 6);
  else if (name == "lost")
    rules = Rules(rules.begin() + 1, rules.end());
  else
    rules = Rules();
  return [rules](Graph& G, Tree& T) {
    leaf_info<Tree> info(G, T);
    int i = 0;
    bool applied = true;
    while(applied && !info.is_path()) {
      applied = false;
      int k = 1;
      for(auto rule : rules) {
        ++k;
        if(rule(G, T, info)) {
          ++i;
          applied = true;
          //std::cerr << ("rule " + std::to_string(k) + "\n");
          //show("step" + std::to_string(i) + ".dot", G, T);
          break;
        }
      }
    }
    return i;
  };
}

template <class Graph>
void run(std::string t, unsigned z, unsigned n, float p, std::string cname, std::string iname) {
  test_suite<Graph> suite(t, z, n, p);
  auto construct = make_construction<Graph>(cname);
  auto improve = make_improvement<Graph,Graph>(iname);
  timing timer;
  unsigned steps;
  double elapsed;
  for(auto G : suite) {
      timer.start();
      auto T = construct(G);
      if(num_vertices(T) != num_vertices(G))
        // G is unconnected
        continue;
      // assert T is a spanning tree
      assert(num_edges(T) == num_vertices(T)-1);
      steps = improve(G, T);
      elapsed = timer.stop();
      // type parameter vertices edges upper construction improvement internal time steps
      std::cout
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
}

int main(int argc, char** argv){
  using std::string;
  using std::vector;

  std::ios_base::sync_with_stdio(0);

  options opt(argc, argv);
  auto z = opt.get<int>("-z", 100);
  auto n = opt.get<int>("-n", 100);
  auto tests = opt.getList<string>("-t", {"gnp", "rgg"});
  auto parameters = opt.getList<float>("-p", {
      0.0001, 0.0005, 0.001, 0.003, 0.005, 0.008, 0.01, 0.03, 0.05, 0.08, 0.1, 0.2, 0.25, 0.3, 0.35,
      //0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99
  });
  auto constructions = opt.getList<string>("-c", {"bfs", "dfs", "rdfs"});
  auto improvements = opt.getList<string>("-i", {"none", "prieto", "lost-light", "lost"});

  //~ vector<float> ps {0.0002, 0.0105, 0.021, 0.0312, 0.0415, 0.0518, 0.062, 0.0725, 0.0827};

  std::cout << "type,parameter,vertices,edges,upper,construction,improvement,internal,time,steps" << std::endl;
  for(auto t : tests)
    for(auto p : parameters)
      for(auto c : constructions)
        for(auto i : improvements)
          run<graph>(t, z, n, p, c, i);

  return 0;
}
