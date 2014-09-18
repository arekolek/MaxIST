// (C) 2014 Arek Olek

#pragma once

#include <algorithm>
#include <vector>

#include <boost/iterator/counting_iterator.hpp>

auto range_iterator = boost::make_counting_iterator<int>;

template<class Graph, class Generator>
void add_spider(Graph& G, unsigned legs, Generator generator) {
  unsigned n = num_vertices(G);
  unsigned cutoff = n / legs;
  std::vector<int> path(range_iterator(0), range_iterator(n));
  std::shuffle(path.begin(), path.end(), generator);
  for(unsigned i = 0; i < n-1; ++i)
    add_edge(path[i % cutoff == 0 ? 0 : i], path[i+1], G);
}

template<class Graph, class Generator>
void add_edges_uniform(Graph& G, double p, Generator generator) {
  std::bernoulli_distribution trial(p);
  unsigned n = num_vertices(G);
  for(unsigned i = 0; i < n; ++i)
    for(unsigned j = i + 1; j < n; ++j)
      if(trial(generator))
        add_edge(i, j, G);
}
