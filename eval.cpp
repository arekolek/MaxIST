// (C) 2014 Arek Olek

#include <iostream>

#include <boost/graph/adjacency_list.hpp>

#include "debug.hpp"
#include "dfs.hpp"
#include "rdfs.hpp"
#include "test_suite.hpp"
#include "options.hpp"

boost::property<boost::edge_color_t,
  boost::default_color_type>            typedef color;
boost::adjacency_list<
  boost::hash_setS, boost::vecS, boost::undirectedS,
  boost::no_property, color>            typedef graph;

template <class Graph>
int eval(Graph const & T) {
  int n = num_vertices(T);
  int internal = 0;
  int upper = n - 2;
  auto vs = vertices(T);
  for(auto vit = vs.first; vit != vs.second; ++vit)
    internal += degree(*vit, T) > 1;

  cout << n-internal << " " << internal << " " << upper
    << " " << n << " " << (upper-internal)/(double)upper << endl;

  return n-internal;
}

int main(int argc, char** argv){
  std::ios_base::sync_with_stdio(0);

  options opt(argc, argv);
  int z = opt.get<int>("-z", 1);
  int n = opt.get<int>("-n", 5);
  float p = opt.get<float>("-p", 1);

  test_suite<graph> suite(z, n, p);

  int sum = 0;

  for(auto G : suite) {
    auto T = dfs_tree(G);

    sum += eval(T);

    T = rdfs_tree(G);

    //show("graph.dot", G, T);
  }

  cout << sum/(double)suite.size() << endl;

  return 0;
}
