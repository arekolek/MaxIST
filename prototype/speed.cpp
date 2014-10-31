// (C) 2014 Arek Olek

#include <iostream>

#include <boost/graph/adjacency_list.hpp>

#include "rdfs.hpp"
#include "range.hpp"

using namespace std;

boost::property<boost::edge_color_t,
  boost::default_color_type>            typedef color;
boost::adjacency_list<
  boost::hash_setS, boost::vecS, boost::undirectedS,
  boost::no_property, color>            typedef graph;


int main() {
  ios_base::sync_with_stdio(0);

  int z, n, m, s, t;
  int internal = 0;

  cin >> z;

  for(int Z = 0; Z < z; ++Z) {
    cin >> n >> m;
    graph g(n);
    for(int i = 0; i < m; ++i) {
      cin >> s >> t;
      add_edge(s, t, g);
    }
    auto T = rdfs_tree(g);
    for(auto v : range(vertices(T)))
      internal += degree(v, T) > 1;
  }
  cout << internal/(double)z << endl;

  return 0;
}
