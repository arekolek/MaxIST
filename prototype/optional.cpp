/*
 * optional.cpp
 *
 *  Created on: Mar 1, 2016
 *      Author: zzzao901
 */

#include <iostream>
#include <vector>
#include <boost/optional.hpp>

using namespace std;
using namespace boost;

optional<size_t> find(int value, const vector<int>& container) {
  for (size_t i = 0; i < container.size(); ++i)
    if (value == container[i])
      return i;
  return boost::none;
}

int main() {
  vector<int> T = {1, 2, 3, 5, 6, 7, 8};
  int value;
  cin >> value;
  if (auto pos = find(value, T)) {
    cout << "found at position " << *pos << endl;
  }
  else {
    cout << "not found\n";
  }
  return 0;
}

