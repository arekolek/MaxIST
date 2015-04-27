// (C) 2014 Arek Olek

#pragma once

#include <algorithm>
#include <vector>
#include <utility>

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

template <class Iterator>
auto shuffled(std::pair<Iterator, Iterator> p) -> std::vector<decltype(*p.first)> {
  std::vector<decltype(*p.first)> R(p.first, p.second);
  std::random_shuffle(R.begin(), R.end());
  return R;
}
