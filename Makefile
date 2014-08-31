
eval.e : eval.cpp $(shell find . -type f -name *.hpp)
	g++ -Wall -O2 -std=c++0x -Igraph -Iutil -Itest -o eval.e eval.cpp
