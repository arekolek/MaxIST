// (C) 2014 Arek Olek

#pragma once

#include <ctime>

class timing {
public:
  double stop() {
    c_end = clock();
    return 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
  }
  void start() {
    c_start = clock();
  }
private:
  clock_t c_start, c_end;
};
