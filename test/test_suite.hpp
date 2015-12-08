// (C) 2014 Arek Olek

#pragma once

#include <chrono>
#include <random>
#include <tuple>

#include <boost/graph/graphml.hpp>

#include "algorithm.hpp"
#include "graph.hpp"

unsigned time_seed() {
  return std::chrono::system_clock::now().time_since_epoch().count();
}

template <class Iterator>
void generate_seeds(Iterator begin, Iterator end) {
  static std::string seed = std::to_string(time_seed());
  std::seed_seq seq(seed.begin(), seed.end());
  seq.generate(begin, end);
}

double pr_within(double x) {
  return x <= 1
    ?  .5*x*x*x*x - 8./3*x*x*x + M_PI*x*x
    : -.5*x*x*x*x - 4*x*x*atan(sqrt(x*x-1)) + 4./3*(2*x*x+1)*sqrt(x*x-1) + (M_PI-2)*x*x + 1./3;
};

template <class Graph>
class test_suite {
public:
  unsigned size() const { return seeds.size() * degrees.size() * sizes.size(); }

  std::string type() const { return t; }

  std::tuple<Graph, double, double> get(unsigned i) const {
    std::default_random_engine generator(seeds[i % seeds.size()]);
    auto Ed = degrees[(i / seeds.size()) % degrees.size()];
    auto n = sizes[i / seeds.size() / degrees.size()];
    auto d = Ed;
    Graph G(n);

    if(t.find("path") != std::string::npos) {
      add_spider(G, 1, generator);
      // the overlap satisfies    y     = a    x     + b
      // full graph satisfies:    0     = a  (n-1)   + b
      // tree graph satisfies: 2(n-1)/n = a 2(n-1)/n + b
      // so we must subtract this:
      d -= 2/(2-double(n)) * d + 1. + double(n)/(n-2);
    }
    if(t.find("rgg") != std::string::npos) {
      auto r = find_argument(d/(n-1), pr_within, 0, sqrt(2));
      Geometric points(n, generator);
      points.add_random_geometric(G, r);
      if(t.find("mst") != std::string::npos) points.add_mst(G);
      return std::make_tuple(G, Ed, r);
    } else {
      auto p = d/(n-1);
      bool mst = t.find("mst") != std::string::npos;
      add_edges_uniform(G, p, generator, mst);
      return std::make_tuple(G, Ed, p);
    }
  }

  template<class Sizes, class Degrees>
  test_suite(std::string t, unsigned z, Sizes ns, Degrees ds) : t(t), sizes(ns), degrees(ds) {
    seeds.resize(z);
    generate_seeds(seeds.begin(), seeds.end());
  }
private:
  std::string t;
  std::vector<unsigned> sizes;
  std::vector<double> degrees;
  std::vector<unsigned> seeds;
};

template <class Graph>
class file_suite {
public:
  std::tuple<Graph, double, double> get(unsigned i) const {
    return std::make_tuple(graphs[i], 0, 0);
  }

  unsigned size() const { return graphs.size(); }

  file_suite(std::string f) : t(f.substr(0, f.find('.'))) {
    std::ifstream file(f);
    if(!file.good()) {
      throw std::invalid_argument("File does not exist: " + f);
    }
    int z, n, m, s, t;
    file >> z;
    while(z--) {
      file >> n >> m;
      Graph G(n);
      for(int i = 0; i < m; ++i) {
        file >> s >> t;
        add_edge(s, t, G);
      }
      graphs.push_back(G);
    }
    file.close();
  }
  std::string type() const {
    return t;
  }
private:
  std::string t;
  std::vector<Graph> graphs;
};

template <class Graph>
class real_suite {
public:
  unsigned size() const { return seeds.size(); }

  std::tuple<Graph, double, double> get(uint i) const {
    std::default_random_engine generator(seeds[i]);
    auto p = shuffled(vertices(G), generator);
    Graph g(num_vertices(G));
    for(auto e : range(edges(G)))
      add_edge(p[source(e, G)], p[target(e, G)], g);
    return std::make_tuple(G, 0, 0);
  }

  real_suite(std::string f, unsigned size) : G(0), t(f.substr(0, f.find('.'))), seeds(size) {
    std::ifstream file(f);
    if(!file.good()) {
      throw std::invalid_argument("File does not exist: " + f);
    }
    boost::dynamic_properties dp;
    read_graphml(file, G, dp);
    generate_seeds(seeds.begin(), seeds.end());
  }
  std::string type() const {
    return t;
  }
private:
  Graph G;
  std::string t;
  std::vector<unsigned> seeds;
};
