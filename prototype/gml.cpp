
#include <fstream>
#include <iostream>
#include <string>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>

#include "range.hpp"
#include "options.hpp"

template<class Graph>
void show(Graph const & g) {
  for(auto v : range(vertices(g))) {
    std::cout << v << ": ";
    for(auto w : range(adjacent_vertices(v, g)))
      std::cout << w << " ";
    std::cout << std::endl;
  }
}

int main(int argc, char* argv[]) {
  options opt(argc, argv);
  auto f = opt.get<std::string>("-f", "");
  std::cout << f << std::endl;

  std::ifstream inFile;
  inFile.open(f, std::ifstream::in);

  typedef boost::adjacency_list<boost::hash_setS, boost::vecS, boost::undirectedS> Graph;
  Graph g;

  boost::dynamic_properties dp;
  read_graphml(inFile, g, dp);

  std::cout << num_vertices(g) << " " << num_edges(g) << std::endl;

  show(g);

  return 0;
}
