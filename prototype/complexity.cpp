// (C) 2014 Arek Olek

#include <iostream>

#include <boost/graph/adjacency_list.hpp>

#include "options.hpp"
#include "test_suite.hpp"
#include "timing.hpp"
#include "bfs.hpp"
#include "lost.hpp"

using namespace std;

typedef boost::adjacency_list<
  boost::hash_setS, boost::vecS, boost::undirectedS
  > graph;
//~ boost::adjacency_matrix<boost::undirectedS> typedef graph;

int main(int argc, char** argv) {
  options opt(argc, argv);
  int z = opt.get<int>("-z", 1);
  int N = opt.get<int>("-n", 500);
  float p = opt.get<float>("-p", 1);

  timing timer;
  lost_light improve;

  for(int n = 100; n <= N; n += 100) {
    test_suite<graph> suite(z, n, p);
    for(auto G : suite) {
      auto T = bfs_tree(G);
      timer.start();
      improve(G, T);
      auto sec = timer.stop();
      cout << "{" << n << "," << sec << "}";
      if(n < N) cout << ",";
      cout << endl;
    }
  }

  return 0;
}
