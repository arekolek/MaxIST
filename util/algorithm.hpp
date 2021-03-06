
#pragma once

#include <algorithm>
#include <cmath>

/**
 * Uses interpolation search to find x such that y-e < f(x) < y+e.
 * The function f must be monotonic.
 **/
template<typename Function>
double find_argument(double y, Function f, double lo, double hi, double e = 1e-9) {
  using std::abs;
  double x;
  while(true) {
    x = lo + (y - f(lo)) * (hi - lo) / (f(hi) - f(lo));

    if(x < lo) return lo;
    if(x > hi) return hi;
    if(abs(f(x)-y) < e) return x;

    if(abs(f(x)-f(lo)) < abs(y-f(lo))) lo = x;
    else hi = x;
  }
}

template <class Container, class Value>
inline unsigned find_index(Container const & container, Value const & value) {
  return std::distance(container.begin(), std::find(container.begin(), container.end(), value));
}

inline bool found(std::string pat, std::string txt) {
  return txt.find(pat) != std::string::npos;
}
