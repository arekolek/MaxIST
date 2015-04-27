// (C) 2014 Arek Olek

#pragma once

#include <ctime>

class timing {
public:
  double stop() {
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &c_end);
    double elapsed = c_end.tv_sec - c_start.tv_sec;
    elapsed += (c_end.tv_nsec - c_start.tv_nsec) / 1000000000.0;
    elapsed *= 1000; // to ms
    return elapsed;
  }
  void start() {
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &c_start);
  }
private:
  timespec c_start, c_end;
};
