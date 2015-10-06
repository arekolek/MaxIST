// (C) 2014 Arek Olek

#pragma once

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <algorithm>
#include <vector>
#include <utility>
#include <random>

template <class Iterator>
class iterator_pair {
public:
  iterator_pair(std::pair<Iterator, Iterator> p_) : p(p_) { }
  Iterator begin() { return p.first; }
  Iterator end() { return p.second; }
private:
  std::pair<Iterator, Iterator> p;
};

template <class Iterator>
iterator_pair<Iterator> range(std::pair<Iterator, Iterator> p) {
  return iterator_pair<Iterator>(p);
}

template <class Iterator>
iterator_pair<Iterator> range(Iterator begin, Iterator end) {
  return iterator_pair<Iterator>(std::make_pair(begin, end));
}

template <class Iterator, class Generator = std::default_random_engine>
auto shuffled(std::pair<Iterator, Iterator> p, Generator generator = Generator()) -> std::vector<decltype(*p.first)> {
  std::vector<decltype(*p.first)> R(p.first, p.second);
  std::shuffle(R.begin(), R.end(), generator);
  return R;
}

template <class IntType = int, class Generator = std::default_random_engine>
IntType random(IntType a = 0, IntType b = std::numeric_limits<IntType>::max()) {
  static Generator generator;
  return std::uniform_int_distribution<IntType>{a, b}(generator);
}

template<class Container>
std::string join(Container const & container, std::string const & delimeter) {
  using boost::algorithm::join;
  using boost::adaptors::transformed;
  using value_type = typename Container::value_type;
  auto tostr = static_cast<std::string(*)(value_type)>(std::to_string);
  return join(container | transformed(tostr), delimeter);
};
