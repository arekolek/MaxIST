// (C) 2014 Arek Olek

#pragma once

#include <algorithm>
#include <chrono>
#include <random>

#include <boost/iterator/counting_iterator.hpp>

auto range_iterator = boost::make_counting_iterator<int>;

unsigned time_seed() {
  return std::chrono::system_clock::now().time_since_epoch().count();
}

template <class Iterator>
void generate_seeds(Iterator begin, Iterator end) {
  std::string seed = std::to_string(time_seed());
  std::seed_seq seq(seed.begin(), seed.end());
  seq.generate(begin, end);
}

template <class Graph>
class test_suite {
public:
  class iterator {
    friend class test_suite;
  public:
    Graph operator *() const;
    const iterator &operator ++() { ++i_; return *this; }
    iterator operator ++(int) { iterator copy(*this); ++i_; return copy; }

    bool operator ==(const iterator &other) const { return i_ == other.i_; }
    bool operator !=(const iterator &other) const { return i_ != other.i_; }

  protected:
    iterator(test_suite const& s, unsigned i) : suite_(s), i_(i) { }

  private:
    test_suite const& suite_;
    unsigned long i_;
  };

  class test_case {
  public:
    unsigned seed() const {
      return s;
    }
    unsigned size() const {
      return n;
    }
    float probability() const {
      return p;
    }
    test_case(unsigned s_, unsigned n_, float p_) : s(s_), n(n_), p(p_) { }
  private:
    unsigned s, n;
    float p;
  };

  iterator begin() const { return iterator(*this, 0); }
  iterator end() const { return iterator(*this, seeds_.size()); }
  test_case test(unsigned i) const { return test_case(seeds_[i], n_, p_); }
  unsigned size() const { return seeds_.size(); }

  test_suite(unsigned z, unsigned n, float p) : n_(n), p_(p) {
    seeds_.resize(z);
    generate_seeds(seeds_.begin(), seeds_.end());
  }
private:
  unsigned n_;
  float p_;
  std::vector<unsigned> seeds_;
};

template <class Graph>
Graph test_suite<Graph>::iterator::operator *() const {
  auto const& test = suite_.test(i_);
  unsigned n = test.size();
  std::default_random_engine generator(test.seed());
  std::bernoulli_distribution trial(test.probability());
  Graph G(n);

  // add random path
  std::vector<int> path(range_iterator(0), range_iterator(n));
  std::shuffle(path.begin(), path.end(), generator);
  for(unsigned i = 0; i < n-1; ++i)
    add_edge(path[i], path[i+1], G);

  // add random edges
  for(unsigned i = 0; i < n; ++i)
    for(unsigned j = i + 1; j < n; ++j)
      if(trial(generator))
        add_edge(i, j, G);

  return G;
}
