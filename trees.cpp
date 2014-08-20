
#include <iostream>
#include <stack>
#include <queue>

#include <boost/graph/adjacency_list.hpp>

#include "debug.hpp"
#include "test_suite.hpp"
#include "range.hpp"

using namespace std;

typedef boost::adjacency_list<
  boost::hash_setS, boost::vecS, boost::undirectedS
  > graph;

typedef pair<int, int> edge;

template <class Graph>
Graph dfs(Graph const & G, int v) {
  Graph T;
  int parent;
  stack<edge> S;
  vector<bool> V(num_vertices(G));
  S.emplace(-1, v);
  while(!S.empty()) {
    tie(parent, v) = S.top();
    S.pop();
    if(!V[v]) {
      V[v] = true;
      if(parent > -1) add_edge(parent, v, T);
      for(auto w : shuffled(adjacent_vertices(v, G))) {
        S.emplace(v, w);
      }
    }
  }
  show("tree-dfs.dot", G, T);
  return T;
}

template <class Graph>
Graph dfs2(Graph const & G, int v) {
  vector<bool> V(num_vertices(G));
  stack<edge> E;
  Graph T;
  while(true) {
    V[v] = true;

    int x = -1, y = -1;
    for(auto w : range(adjacent_vertices(v, G)))
      if(!V[w]) {
        tie(x, y) = tie(v, w);
        E.emplace(v, w);
      }

    if(y == -1) {
      while(!E.empty() && V[E.top().second]) E.pop();
      if(E.empty()) break;
      tie(x, y) = E.top();
    } else {
      E.pop();
    }

    add_edge(x, y, T);

    v = y;
  }
  return T;
}

template <class Graph>
Graph dfs_fifo(Graph const & G, int v) {
  int rank;
  vector<int> V(num_vertices(G), rank = 0);
  queue<edge> E;
  Graph T;
  while(true) {
    V[v] = ++rank;

    int x, y = -1;
    for(auto w : range(adjacent_vertices(v, G)))
      if(V[w] == 0) {
        tie(x, y) = tie(v, w);
        break;
      }

    if(y == -1) {
      while(!E.empty() && V[E.front().second] != 0) E.pop();
      if(E.empty()) break;
      tie(x, y) = E.front();
    }

    for(int w : range(adjacent_vertices(x, G)))
      if(w != y && V[w] == 0) {
        E.emplace(x, w);
      }

    add_edge(x, y, T);

    v = y;
  }
  show("tree-dfs-fifo.dot", G, T);
  return T;
}

template <class Graph>
Graph dfs_fifo2(Graph const & G, int v) {
  vector<bool> V(num_vertices(G));
  queue<edge> E;
  Graph T;
  while(true) {
    V[v] = true;

    int x = -1, y = -1;
    for(auto w : range(adjacent_vertices(v, G)))
      if(!V[w]) {
        if(y == -1) tie(x, y) = tie(v, w);
        else E.emplace(v, w);
      }

    if(y == -1) {
      while(!E.empty() && V[E.front().second]) E.pop();
      if(E.empty()) break;
      tie(x, y) = E.front();
    }

    add_edge(x, y, T);

    v = y;
  }
  return T;
}

template <class Graph>
Graph bfs(Graph const & G, int v) {
  Graph T;
  int parent;
  queue<edge> Q;
  vector<bool> V(num_vertices(G));
  Q.emplace(-1, v);
  V[v] = true;
  while(!Q.empty()) {
    tie(parent, v) = Q.front();
    Q.pop();
    for(auto w : range(adjacent_vertices(v, G))) {
      if(!V[w]) {
        V[w] = true;
        add_edge(v, w, T);
        Q.emplace(v, w);
      }
    }
  }
  show("tree-bfs.dot", G, T);
  return T;
}

bool operator==(graph const & A, graph const & B) {
  if(num_vertices(A) != num_vertices(B))
    return false;
  for(auto e : range(edges(A)))
    if(!boost::edge(source(e, A), target(e, A), B).second)
      return false;
  return true;
}

bool operator!=(graph const & A, graph const & B) {
  return !(A == B);
}

int main() {
  std::srand ( unsigned ( std::time(0) ) );

  ios_base::sync_with_stdio(0);

  test_suite<graph> suite(1, 20, 0.2);

  for(auto const& G : suite) {
    auto dfst = dfs(G, 0);
    auto dfst2 = dfs(G, 0);
    auto fdfst = dfs_fifo(G, 0);
    auto fdfst2 = dfs_fifo2(G, 0);
    auto bfst = bfs(G, 0);
    cout << (dfst == dfst2 ? "tak" : "nie") << endl;
    if(dfst != dfst2) {
      show("rozny1.dot", G, dfst);
      show("rozny2.dot", G, dfst2);
    }
  }

  return 0;
}