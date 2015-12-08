// (C) 2014 Arek Olek

#pragma once

#include <chrono>
#include <random>

#include <boost/graph/graphml.hpp>

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

template <class Graph>
class test_suite {
public:
  unsigned size() const { return seeds.size(); }

  float parameter() const { return p; }
  std::string type() const { return t; }

  Graph get(unsigned i) const {
    Graph G(n);
    std::default_random_engine generator(seeds[i]);
    if(t.find("path") != std::string::npos) add_spider(G, 1, generator);
    if(t.find("rgg") != std::string::npos) {
      Geometric points(n, generator);
      points.add_random_geometric(G, p);
      if(t.find("mst") != std::string::npos) points.add_mst(G);
    }
    if(t.find("gnp") != std::string::npos) {
      bool mst = t.find("mst") != std::string::npos;
      add_edges_uniform(G, p, generator, mst);
    }
    return G;
  }

  test_suite(std::string t, unsigned z, unsigned n, float p) : t(t), n(n), p(p) {
    seeds.resize(z);
    generate_seeds(seeds.begin(), seeds.end());
  }
private:
  std::string t;
  unsigned n;
  float p;
  std::vector<unsigned> seeds;
};

template <class Graph>
class file_suite {
public:
  class iterator {
    friend class file_suite;
  public:
    Graph operator *() const {
      int n, m, s, t;
      suite.file >> n >> m;
      Graph G(n);
      for(int i = 0; i < m; ++i) {
        suite.file >> s >> t;
        add_edge(s, t, G);
      }
      return G;
    }
    const iterator &operator ++() { ++i; return *this; }
    iterator operator ++(int) { iterator copy(*this); ++i; return copy; }

    bool operator ==(const iterator &other) const { return i == other.i; }
    bool operator !=(const iterator &other) const { return i != other.i; }

  protected:
    iterator(file_suite const& s, unsigned i) : suite(s), i(i) {}

  private:
    file_suite const& suite;
    unsigned long i;
  };

  iterator begin() const { return iterator(*this, 0); }
  iterator end() const { return iterator(*this, size_); }
  unsigned size() const { return size_; }

  Graph get(uint i) const {
    return *iterator(*this, i);
  }

  file_suite(std::string f) : t(f.substr(0, f.find('.'))), file(f), size_(0) {
    if(!file.good()) {
      throw std::invalid_argument("File does not exist: " + f);
    }
    file >> size_;
  }
  std::string type() const {
    return t;
  }
  int parameter() const {
    return 0;
  }
private:
  std::string t;
  mutable std::ifstream file;
  unsigned size_;
};

template <class Graph>
class real_suite {
public:
  unsigned size() const { return seeds.size(); }

  Graph get(uint i) const {
    std::default_random_engine generator(seeds[i]);
    auto p = shuffled(vertices(G), generator);
    Graph g(num_vertices(G));
    for(auto e : range(edges(G)))
      add_edge(p[source(e, G)], p[target(e, G)], g);
    return g;
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
  int parameter() const {
    return 0;
  }
private:
  Graph G;
  std::string t;
  std::vector<unsigned> seeds;
};
