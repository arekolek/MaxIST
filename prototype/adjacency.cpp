
#include <fstream>
#include <iostream>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/random.hpp>

#include "range.hpp"

using namespace std;

template<class Graph, class Generator>
Graph random_graph(int n, double p, Generator generator) {
  Graph g(n);
  bernoulli_distribution distribution(p);
  auto trial = bind(distribution, generator);
  for(int i = 0; i < n; ++i)
    for(int j = i + 1; j < n; ++j)
      if(trial()) add_edge(i, j, g);
  return g;
}

class Timing {
  timespec c_start, c_end;
  clockid_t c_id;
  vector<time_t> times;
public:
  Timing(clockid_t id) : c_id(id) {}
  void start() {
    clock_gettime(c_id, &c_start);
  }
  void stop() {
    clock_gettime(c_id, &c_end);
  }
  void tick() {
    stop();
    times.push_back(elapsed());
    start();
  }
  time_t begin() const {
    return c_start.tv_sec * 1e9 + c_start.tv_nsec;
  }
  time_t end() const {
    return c_end.tv_sec * 1e9 + c_end.tv_nsec;
  }
  time_t elapsed() const {
    return end() - begin();
  }
  const vector<time_t>& elems() {
    return times;
  }
  void clear() {
    times.clear();
  }
};

template <class Graph>
void test(string name) {
  default_random_engine gen;
  ofstream file(name);
  Timing t(CLOCK_PROCESS_CPUTIME_ID);
  for(int n = 100; n <= 1000; n += 100) {
    cerr << n << endl;
    for(int j = 0; j < 100; ++j) {
      auto g = random_graph<Graph>(n, 0.1, gen);

      int foo = 0;
      t.start();
      for(auto v : range(vertices(g))) ++foo;
      t.stop();
      auto t1 = t.elapsed();

      t.start();
      for(auto e : range(edges(g))) --foo;
      t.stop();
      auto t2 = t.elapsed();

      for(int i = 0; i < 10; ++i) {
        auto v = random_vertex(g, gen);
        auto w = random_vertex(g, gen);

        t.start();
        auto e = edge(v, w, g);
        t.tick();
        auto d = out_degree(v, g);
        t.tick();
        for(auto x : range(adjacent_vertices(v, g))) --d;
        t.tick();

        if(e.second) remove_edge(v, w, g);

        t.start();
        add_edge(v, w, g);
        t.tick();
        remove_edge(v, w, g);
        t.tick();

        if(e.second) add_edge(v, w, g);

        t.start();
        source(e.first, g);
        t.tick();
        target(e.first, g);
        t.tick();

        assert(d == 0 && foo < 0);

        file << n << ' ' << t1 << ' ' << t2 << ' ';
        for(auto x : t.elems()) file << x << ' ';
        file << endl;
        t.clear();
      }
    }
  }
  file.close();
}

int main() {
  using namespace boost;
  test<adjacency_list<vecS, vecS, undirectedS>>("vec.txt");
  test<adjacency_list<listS, vecS, undirectedS>>("list.txt");
  test<adjacency_list<slistS, vecS, undirectedS>>("slist.txt");
  test<adjacency_list<setS, vecS, undirectedS>>("set.txt");
  test<adjacency_list<multisetS, vecS, undirectedS>>("multiset.txt");
  test<adjacency_list<hash_setS, vecS, undirectedS>>("hashset.txt");
  test<adjacency_matrix<undirectedS>>("matrix.txt");
  return 0;
}
