// (C) 2014 Arek Olek

#include <functional>
#include <iostream>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "bfs.hpp"
#include "dfs.hpp"
#include "rdfs.hpp"
#include "lost.hpp"

#include "debug.hpp"
#include "options.hpp"
#include "timing.hpp"
#include "range.hpp"

#include "test_suite.hpp"

using namespace std;

boost::property<boost::edge_color_t,
  boost::default_color_type>            typedef color;
boost::adjacency_list<
  boost::hash_setS, boost::vecS, boost::undirectedS,
  boost::no_property, color>            typedef graph;

template <class Graph>
int eval(Graph const & T) {
  int internal = 0;
  for(auto v : range(vertices(T)))
    internal += degree(v, T) > 1;
  return internal;
}

template <class Graph>
float average_degree(Graph const & G) {
  float sum = 0;
  for(auto v : range(vertices(G)))
    sum += degree(v, G);
  return sum / num_vertices(G);
}

function<graph(graph&)> typedef solution;

class dummy {
public:
  template<class Graph>
  int operator()(Graph& G, Graph& T) {
    return 0;
  }
};

template <class Improvement>
void run(int z, int n, vector<float> ps, vector<string> name,
  vector<solution> algo, Improvement improve, string iname) {
  timing timer;

  cout << endl << iname << endl;

  cout << "\t\t";
  for(auto n : name)
    cout << n << "\t\t\t\t\t";
  cout << "\np\tE[deg]";
  for(auto n : name)
    cout << "\tmin\tE[|I|]\tmax\ttime\tsteps";
  cout << endl;

  for(auto p : ps) {
    test_suite<graph> suite(z, n, p);

    double degree = 0;
    vector<double>
      steps(algo.size(), 0),
      quality(algo.size(), 0),
      time(algo.size(), 0);
    vector<int>
      minimum(algo.size(), n),
      maximum(algo.size(), 0);

    for(auto G : suite) {
      degree += average_degree(G);

      for(unsigned i = 0; i < algo.size(); ++i) {
        timer.start();
        auto T = algo[i](G);
        steps[i] += improve(G, T);
        time[i] += timer.stop();
        int internals = eval(T);
        quality[i] += internals;
        minimum[i] = min(minimum[i], internals);
        maximum[i] = max(maximum[i], internals);
        //show("graph" + to_string(i) + ".dot", G, T);
      }
    }

    int count = suite.size();
    cout << p << '\t' << degree / count;

    for(unsigned i = 0; i < algo.size(); ++i)
      cout
        << '\t' << minimum[i]
        << '\t' << quality[i] / count
        << '\t' << maximum[i]
        << '\t' << time[i] / count
        << '\t' << steps[i] / count
        ;

    cout << endl;
  }
}

int main(int argc, char** argv){
  ios_base::sync_with_stdio(0);

  options opt(argc, argv);
  int z = opt.get<int>("-z", 1);
  int n = opt.get<int>("-n", 5);
  float p = opt.get<float>("-p", -1);
  string a = opt.get<string>("-a");
  vector<float> ps {0.0001, 0.0005, 0.001, 0.003, 0.005, 0.008,
    0.01, 0.03, 0.05, 0.08,
    0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
    0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99};
  if(p >= 0 && p <= 1) {
    ps.clear();
    ps.push_back(p);
  }

  vector<string> name {"bfs", "dfs", "rdfs"};
  vector<solution> algo {bfs_tree<graph>, dfs_tree<graph>, rdfs_tree<graph>};

  run(z, n, ps, name, algo, dummy(), "no improvement");

  run(z, n, ps, name, algo, prieto(), "prieto");

  run(z, n, ps, name, algo, lost_light(), "lost-light");

  return 0;
}
