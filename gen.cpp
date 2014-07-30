
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

using namespace std;

using boost::adjacency_list;
using boost::graph_traits;

adjacency_list<
  boost::hash_setS, boost::vecS,
  boost::undirectedS>                   typedef Graph;

graph_traits<Graph>::vertex_descriptor  typedef Vertex;
graph_traits<Graph>::vertex_iterator    typedef VertexIterator;
graph_traits<Graph>::edge_descriptor    typedef Edge;
graph_traits<Graph>::out_edge_iterator  typedef EdgeIterator;
graph_traits<Graph>::edge_iterator      typedef AllEdgesIterator;

int getiopt(int argc, char** argv, string opt, int def = 0){
  for(int i = 1; i < argc-1; ++i) if(strcmp(argv[i], opt.c_str()) == 0) return atoi(argv[i+1]);
  return def;
}
float getfopt(int argc, char** argv, string opt, float def = 0){
  for(int i = 1; i < argc-1; ++i) if(strcmp(argv[i], opt.c_str()) == 0) return atof(argv[i+1]);
  return def;
}

struct graph_writer {
  void operator()(std::ostream& out) const {
    out << "node [shape=point color=red]" << endl;
    out << "edge [style=dashed]" << endl;
  }
};

template<class G, class GPW>
void graphviz(string file, G g, GPW gpw) {
  ofstream f(file);
  boost::write_graphviz(f, g,
    boost::default_writer(), boost::default_writer(), gpw);
  f.close();
}

int main(int argc, char** argv) {
  int z = getiopt(argc, argv, "-z", 1);
  int n = getiopt(argc, argv, "-n", 5);
  float p = getfopt(argc, argv, "-p", 1);

  cout << z << endl;
  while(z--){

    Graph G(n);

    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    bernoulli_distribution distribution(p);

    for(int i = 0; i < n; ++i)
      for(int j = i + 1; j < n; ++j)
        if(distribution(generator))
          boost::add_edge(i, j, G);

    vector<int> path(n);
    for(int i = 0; i < n; ++i) path[i] = i;
    shuffle(path.begin(), path.end(), generator);
    for(int i = 0; i < n-1; ++i)
      boost::add_edge(path[i], path[i+1], G);

    cout << n << " " << boost::num_edges(G) << endl;
    auto es = edges(G);
    for(auto eit = es.first; eit != es.second; ++eit)
      cout << source(*eit, G) << " " << target(*eit, G) << endl;

    //graphviz("graph.dot", G, graph_writer());

  }

  return 0;
}
