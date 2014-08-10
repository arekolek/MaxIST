// (C) 2014 Arek Olek

#pragma once

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
