// (C) 2014 Arek Olek

#pragma once

#include <chrono>
#include <random>

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
  class iterator {
    friend class test_suite;
  public:
    Graph operator *() const {
      Graph G(suite.num_vertices());
      std::default_random_engine generator(suite.seed(i));
      add_spider(G, 1, generator);
      if(suite.type() == "rgg")
        add_random_geometric(G, suite.parameter(), generator);
      else
        add_edges_uniform(G, suite.parameter(), generator);
      return G;
    }
    const iterator &operator ++() { ++i; return *this; }
    iterator operator ++(int) { iterator copy(*this); ++i; return copy; }

    bool operator ==(const iterator &other) const { return i == other.i; }
    bool operator !=(const iterator &other) const { return i != other.i; }

  protected:
    iterator(test_suite const& s, unsigned i) : suite(s), i(i) { }

  private:
    test_suite const& suite;
    unsigned long i;
  };

  iterator begin() const { return iterator(*this, 0); }
  iterator end() const { return iterator(*this, size()); }
  unsigned size() const { return seeds.size(); }

  unsigned seed(unsigned i) const { return seeds[i]; }
  unsigned num_vertices() const { return n; }
  float parameter() const { return p; }
  std::string type() const { return t; }

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
  iterator end() const { return iterator(*this, size); }

  file_suite(std::string f) : t(f.substr(0, f.find('.'))), file(f), size(0) {
    if(!file.good()) {
      throw std::invalid_argument("File does not exist: " + f);
    }
    file >> size;
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
  unsigned size;
};
