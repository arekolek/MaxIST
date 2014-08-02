
#include <iostream>

#include <boost/graph/adjacency_list.hpp>

#include "options.hpp"
#include "test_suite.hpp"

using namespace std;

typedef boost::adjacency_list<
  boost::hash_setS, boost::vecS, boost::undirectedS
  > graph;

int main(int argc, char** argv) {
  options opt(argc, argv);
  int z = opt.get<int>("-z", 1);
  int n = opt.get<int>("-n", 5);
  float p = opt.get<float>("-p", 1);

  test_suite<graph> suite(z, n, p);

  cout << suite.size() << endl;

  for(auto const& G : suite) {

    cout << num_vertices(G) << " " << num_edges(G) << endl;
    auto es = edges(G);
    for(auto eit = es.first; eit != es.second; ++eit)
      cout << source(*eit, G) << " " << target(*eit, G) << endl;

  }

  return 0;
}
