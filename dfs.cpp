// (C) 2014 Arek Olek

#include <iostream>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/undirected_dfs.hpp>

using std::cin;
using std::cout;
using std::endl;
using std::vector;

using boost::adjacency_list;
using boost::edge_color;
using boost::graph_traits;
using boost::make_dfs_visitor;
using boost::num_vertices;
using boost::on_tree_edge;
using boost::record_predecessors;
using boost::write_graphviz;
using boost::vertices;
using boost::visitor;
using boost::undirected_dfs;

boost::property<boost::edge_color_t,
  boost::default_color_type>            typedef Color;
adjacency_list<
  boost::hash_setS, boost::vecS, boost::undirectedS,
  boost::no_property, Color>            typedef Graph;
boost::property_map<
  Graph, boost::edge_color_t>::type     typedef ColorMap;

graph_traits<Graph>::vertex_descriptor  typedef Vertex;
graph_traits<Graph>::vertex_iterator    typedef VertexIterator;
graph_traits<Graph>::edge_descriptor    typedef Edge;
graph_traits<Graph>::out_edge_iterator  typedef EdgeIterator;
graph_traits<Graph>::edge_iterator      typedef AllEdgesIterator;

template <class ParentDecorator>
struct print_parent {
  print_parent(const ParentDecorator& p_) : p(p_) { }
  template <class Vertex>
  void operator()(const Vertex& v) const {
    std::cout << "parent[" << v << "] = " <<  p[v]  << std::endl;
  }
  ParentDecorator p;
};


template <class NewGraph, class Tag>
struct graph_copier
  : public boost::base_visitor<graph_copier<NewGraph, Tag> >
{
  typedef Tag event_filter;

  graph_copier(NewGraph& graph) : new_g(graph) { }

  template <class Edge, class Graph>
  void operator()(Edge e, Graph& g) {
    boost::add_edge(boost::source(e, g), boost::target(e, g), new_g);
  }
private:
  NewGraph& new_g;
};

template <class NewGraph, class Tag>
inline graph_copier<NewGraph, Tag>
copy_graph(NewGraph& g, Tag) {
  return graph_copier<NewGraph, Tag>(g);
}

struct red_points_dashed {
  void operator()(std::ostream& out) const {
    out << "node [shape=point color=brown]" << endl;
    out << "edge [style=dashed color=burlywood1]" << endl;
  }
};

template <class G>
class node_writer {
public:
  node_writer(G const& t) : tree(t) {}
  template <class V>
  void operator()(std::ostream& out, const V& v) const {
    if(boost::degree(v, tree) == 1)
      out << "[style=solid color=green]";
    if(boost::degree(v, tree) > 2)
      out << "[style=solid color=red]";
  }
private:
  G const& tree;
};

template <class G>
node_writer<G> make_node_writer(G t) {
  return node_writer<G>(t);
}

template <class G>
class edge_writer {
public:
  edge_writer(G const& g, G const& t) : graph(g), tree(t) {}
  template <class E>
  void operator()(std::ostream& out, const E& e) const {
    auto vs = adjacent_vertices(source(e, graph), tree);
    for(auto vit = vs.first; vit != vs.second; ++vit)
      if(*vit == target(e, graph)) {
        out << "[style=solid color=burlywood]";
        break;
      }
  }
private:
  G const& graph;
  G const& tree;
};

template <class G>
edge_writer<G> make_edge_writer(G g, G t) {
  return edge_writer<G>(g, t);
}

template<class G>
void show(std::string file, G g, G t) {
  std::ofstream f(file);
  boost::write_graphviz(f, g,
    make_node_writer(t), make_edge_writer(g, t), red_points_dashed());
  f.close();
}

int main(int argc, char** argv){
  std::ios_base::sync_with_stdio(0);

  int Z, sum = 0;
  cin >> Z;

  for(int z = 0; z < Z; ++z) {
    int n, m;
    cin >> n >> m;

    Graph G(n);
    for(int i = 0; i < m; ++i) {
      Vertex a, b;
      cin >> a >> b;
      add_edge(a, b, G);
    }

    vector<Vertex> p(num_vertices(G));
    for(unsigned v = 0; v < p.size(); ++v) p[v] = v;

    Graph T;

    undirected_dfs(G,
      visitor(make_dfs_visitor(copy_graph(T, on_tree_edge())))
        .edge_color_map(get(edge_color, G)));

    int internal = 0;
    int upper = n - 2;
    auto vs = vertices(T);
    for(auto vit = vs.first; vit != vs.second; ++vit)
      internal += boost::degree(*vit, T) > 1;

    cout << n-internal << " " << internal << " " << upper
      << " " << n << " " << (upper-internal)/(double)upper << endl;

    sum += n-internal;

    show("graph.dot", G, T);

  }

  cout << sum/(double)Z << endl;

  return 0;
}

