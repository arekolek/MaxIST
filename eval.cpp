// (C) 2014 Arek Olek

#include <functional>
#include <iostream>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "debug.hpp"
#include "dfs.hpp"
#include "rdfs.hpp"
#include "test_suite.hpp"
#include "options.hpp"
#include "timing.hpp"

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

  cerr << n-internal << " " << internal << " " << upper
    << " " << n << " " << (upper-internal)/(double)upper << endl;

  return n-internal;
}

template <class Graph>
float average_degree(Graph const & G) {
  auto vs = vertices(G);
  float sum = 0;
  for(auto vit = vs.first; vit != vs.second; ++vit)
    sum += degree(*vit, G);
  return sum / num_vertices(G);
}

function<graph(graph&)> typedef solution;

int main(int argc, char** argv){
  ios_base::sync_with_stdio(0);

  options opt(argc, argv);
  int z = opt.get<int>("-z", 1);
  int n = opt.get<int>("-n", 5);
  //float p = opt.get<float>("-p", 1);
  auto ps = {0.0001, 0.0005, 0.001, 0.003, 0.005, 0.008,
    0.01, 0.03, 0.05, 0.08,
    0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
    0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99};

  vector<solution> algo {dfs_tree<graph>, rdfs_tree<graph>};

  timing timer;

  for(auto p : ps) {
    test_suite<graph> suite(z, n, p);

    double degree = 0;
    vector<double> quality(algo.size(), 0), time(algo.size(), 0);

    for(auto G : suite) {
      degree += average_degree(G);

      for(unsigned i = 0; i < algo.size(); ++i) {
        timer.start();
        auto T = algo[i](G);
        time[i] += timer.stop();
        quality[i] += eval(T);
        cerr << endl;
      }

      //show("graph.dot", G, T);
    }

    int count = suite.size();
    cout << p << '\t' << degree / count << '\t';

    for(unsigned i = 0; i < algo.size(); ++i)
      cout
        << quality[i] / count << '\t'
        << time[i] / count << '\t';

    cout << endl;
  }


  return 0;
}
