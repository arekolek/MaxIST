/*
 * for_loop.cpp
 *
 *  Created on: Mar 1, 2016
 *      Author: zzzao901
 */

#include <iostream>
#include <chrono>
#include <thread>

using namespace std;

int foo() {
  this_thread::sleep_for(std::chrono::milliseconds(1000));
  return 10;
}

int main(int argc, char **argv) {
  for(int i = 0; i < foo(); ++i) {
    cout << "Elo\n";
  }
  return 0;
}
