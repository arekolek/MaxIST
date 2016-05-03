
#include <random>

#include <boost/graph/adjacency_list.hpp>

#include "graph.hpp"
#include "debug.hpp"

boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> typedef avec;

int main(int argc, char **argv) {
  std::random_device rd;
  std::default_random_engine generator(rd());
  avec g(20);
  add_spider(g, 3, generator);
  show("path.dot", g);
  return 0;
}
