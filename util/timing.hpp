// (C) 2014 Arek Olek

#pragma once

#include <ctime>
#include <string>

class timing {
  timespec c_start, c_end;
  clockid_t c_id;
public:
  timing(clockid_t id = CLOCK_THREAD_CPUTIME_ID) : c_id(id) {}
  double stop() {
    clock_gettime(c_id, &c_end);
    double elapsed = (c_end.tv_sec - c_start.tv_sec) * 1000.0;
    elapsed += (c_end.tv_nsec - c_start.tv_nsec) / 1000000.0;
    return elapsed;
  }
  void start() {
    clock_gettime(c_id, &c_start);
  }
};

std::string readable(int t) {
  if(t < 1000) return std::to_string(t) + " ms";
  t /= 1000;

  std::string result = std::to_string(t % 60) + " s";
  if(t < 60) return result;
  t /= 60;

  result = std::to_string(t % 60) + " m " + result;
  if(t < 60) return result;
  t /= 60;

  return std::to_string(t) + " h " + result;
}
