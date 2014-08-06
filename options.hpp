// (C) 2014 Arek Olek

#pragma once

#include <cstring>
#include <string>

namespace opt_util {
  template <class Result>
  Result convert(char* cstr) {
    throw "stub";
  }
  template <>
  int convert<int>(char* cstr) {
    return atoi(cstr);
  }
  template <>
  float convert<float>(char* cstr) {
    return atof(cstr);
  }
  template <>
  std::string convert<std::string>(char* cstr) {
    return std::string(cstr);
  }
}

class options {
public:
  template <class Result>
  Result get(std::string name, Result def = Result()) {
    for(int i = 1; i < argc-1; ++i)
      if(strcmp(argv[i], name.c_str()) == 0)
        return opt_util::convert<Result>(argv[i+1]);
    return def;
  }
  options(int c, char** v) : argc(c), argv(v) { }
private:
  int argc;
  char** argv;
};
