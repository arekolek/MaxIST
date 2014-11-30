// (C) 2014 Arek Olek

#pragma once

#include <cstring>
#include <string>
#include <sstream>
#include <vector>

namespace opt_util {
  template <class Result>
  Result convert(const char* cstr) {
    throw "stub";
  }
  template <>
  int convert<int>(const char* cstr) {
    return atoi(cstr);
  }
  template <>
  float convert<float>(const char* cstr) {
    return atof(cstr);
  }
  template <>
  std::string convert<std::string>(const char* cstr) {
    return std::string(cstr);
  }
}

class options {
public:
  unsigned find(std::string const& name) const {
    for(int i = 1; i < argc-1; ++i)
      if(strcmp(argv[i], name.c_str()) == 0)
        return i + 1;
    return 0;
  }
  template <class Result = std::string>
  Result get(std::string name, Result def = Result()) {
    return find(name) ? opt_util::convert<Result>(argv[find(name)]) : def;
  }
  template <class Result>
  std::vector<Result> getList(std::string name, std::initializer_list<Result> def = {}) {
    std::vector<Result> result;
    std::stringstream text(get(name));
    for(std::string token; std::getline(text, token, ',');)
      result.push_back(opt_util::convert<Result>(token.c_str()));
    return find(name) ? result : def;
  }
  options(int c, char** v) : argc(c), argv(v) { }
private:
  int argc;
  char** argv;
};
