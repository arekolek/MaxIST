
#include <functional>
#include <iostream>
#include <random>

using namespace std;

template<class Generator>
double foo(Generator& generator) {
  uniform_real_distribution<> distribution(0, 1);
  auto trial = bind(distribution, ref(generator));
  return trial();
}

int main(int argc, char **argv) {
  random_device rd;
  default_random_engine generator(rd());
  cout << foo(generator) << endl;
  cout << foo(generator) << endl;
  return 0;
}
