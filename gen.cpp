
#include <cstring>
#include <iostream>

#include <boost/graph/adjacency_list.hpp>

#include "test_suite.hpp"

using namespace std;

typedef boost::adjacency_list<
  boost::hash_setS, boost::vecS, boost::undirectedS
  > graph;

int getiopt(int argc, char** argv, string opt, int def = 0){
  for(int i = 1; i < argc-1; ++i) if(strcmp(argv[i], opt.c_str()) == 0) return atoi(argv[i+1]);
  return def;
}
float getfopt(int argc, char** argv, string opt, float def = 0){
  for(int i = 1; i < argc-1; ++i) if(strcmp(argv[i], opt.c_str()) == 0) return atof(argv[i+1]);
  return def;
}

int main(int argc, char** argv) {
  int z = getiopt(argc, argv, "-z", 1);
  int n = getiopt(argc, argv, "-n", 5);
  float p = getfopt(argc, argv, "-p", 1);

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
