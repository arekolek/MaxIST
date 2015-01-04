// (C) 2014 Arek Olek

#include <iostream>

#include <boost/graph/adjacency_list.hpp>

#include "bfs.hpp"
#include "lost.hpp"

#include "test_suite.hpp"

#include "options.hpp"
#include "timing.hpp"

using namespace std;

typedef boost::adjacency_list<
  boost::hash_setS, boost::vecS, boost::undirectedS
  > list_graph;
boost::adjacency_matrix<boost::undirectedS> typedef graph;

int main(int argc, char** argv) {
  options opt(argc, argv);
  int z = opt.get<int>("-z", 1);
  int N = opt.get<int>("-n", 500);
  float p = opt.get<float>("-p", 1);

  timing timer;
  auto improve = make_improvement<graph, graph>("prieto");

  for(int n = 100; n <= N; n += 100) {
    test_suite<graph> suite("gnp", z, n, p);
    double total = 0;
    for(auto G : suite) {
      auto T = bfs_tree(G);
      //list_graph E;
      //copy_edges(T, E);
      timer.start();
      improve(G, T);
      total += timer.stop();
    }
    cout << "{" << n << "," << total/suite.size() << "}";
    if(n < N) cout << ",";
    cout << endl;
  }

  return 0;
}
