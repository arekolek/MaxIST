// (C) 2014 Arek Olek

#include <random>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "debug.hpp"
#include "options.hpp"
#include "range.hpp"

template <class Vertex, class Graph>
Vertex random_neighbor(Vertex const & v, Graph const & G) {
  unsigned target = random<unsigned>(0, out_degree(v, G)-1);
  for(auto w : range(adjacent_vertices(v, G)))
    if(target-- == 0) return w;
  assert(false);
}

template <class Graph>
unsigned random_walk(Graph const & G, unsigned v, unsigned t) {
  unsigned steps = 0;
  while(v != t) {
    v = random_neighbor(v, G);
    ++steps;
  }
  return steps;
}

template<class Graph>
void barbell(Graph& G) {
  auto n = num_vertices(G);
  auto m = n / 3;
  for(int i = 0; i < m; ++i) {
    for(int j = i+1; j < m; ++j) {
      add_edge(i, j, G);
      add_edge(i+2*m, j+2*m, G);
    }
    add_edge(i+m, i+m+1, G);
  }
  add_edge(m-1, m, G);
}

boost::adjacency_list<boost::hash_setS, boost::vecS, boost::undirectedS> typedef alist;

int main(int argc, char** argv) {
  options opt(argc, argv);
  const auto z = opt.get<int>("-z", 100);
  auto n = 3*opt.get<int>("-m", 5);
  alist g(n);
  barbell(g);
  std::vector<int> degrees;
  for(auto v : range(vertices(g))) degrees.push_back(degree(v, g));
  std::default_random_engine gen;
  std::discrete_distribution<> pi(degrees.begin(), degrees.end());
  for(int i = 0; i < z; ++i)
    std::cout << n << '\t' << random_walk(g, pi(gen), pi(gen)) << std::endl;
  return 0;
}
