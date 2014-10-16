
includes = -Igraph -Iutil -Itest
CXXFLAGS = -Wall -O2 -std=c++0x $(includes)

eval.e : eval.cpp $(shell find . -type f -name *.hpp)
	$(CXX) $(CXXFLAGS) -o $@ $<
