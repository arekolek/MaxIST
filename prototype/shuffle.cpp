
#include <algorithm>
#include <iostream>
#include <vector>
#include <random>

using namespace std;

template <class Container, class Generator = std::default_random_engine>
Container test(Container const & container) {
  Container R(container);
  std::random_shuffle(R.begin(), R.end());
  return R;
}

int main() {
  vector<int> a{1,2,3,4,5,6,7};
  for(auto i : a) cout << i << ' ';
  cout << endl;
  for(auto i : test(a)) cout << i << ' ';
  cout << endl;
  for(auto i : test(a)) cout << i << ' ';
  cout << endl;
  for(auto i : a) cout << i << ' ';
  cout << endl;
}
