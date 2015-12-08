
#include <iostream>

#include "algorithm.hpp"

using namespace std;

int main() {
  auto f = [](double x) {
    return x<=1
      ?  .5*x*x*x*x - 8./3*x*x*x + M_PI*x*x
      : -.5*x*x*x*x - 4*x*x*atan(sqrt(x*x-1)) + 4./3*(2*x*x+1)*sqrt(x*x-1) + (M_PI-2)*x*x + 1./3; };
  double y = 0.99;
  auto x = find_argument(y, f, 0, 1.4);
  cout << "f(" << x << ") = " << f(x) << " ≈ " << y << endl;
  return 0;
}
