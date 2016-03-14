// (C) 2014 Arek Olek

#pragma once

#include <boost/algorithm/string/join.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <algorithm>
#include <vector>
#include <utility>
#include <random>

template <typename ... Args>
auto range(Args&& ... args) -> decltype(boost::make_iterator_range(std::forward<Args>(args)...)) {
   return boost::make_iterator_range(std::forward<Args>(args)...);
}

template <typename ... Args>
auto filter(Args&& ... args) {
   return boost::make_filter_iterator(std::forward<Args>(args)...);
}

auto range_iterator(int a, int b) {
  return std::make_pair(boost::make_counting_iterator<int>(a),
                        boost::make_counting_iterator<int>(b));
}

template <class Iterator, class Generator = std::default_random_engine>
auto shuffled(std::pair<Iterator, Iterator> p, Generator generator = Generator()) {
  std::vector<typename std::iterator_traits<Iterator>::value_type> R(p.first, p.second);
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

template<class Container>
int sum(Container const & container) {
  return std::accumulate(container.begin(), container.end(), 0);
};
