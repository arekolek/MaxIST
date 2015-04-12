// (C) 2014 Arek Olek

#pragma once

template <class Graph>
Graph five_three_tree(Graph const & G) {
  const int n = num_vertices(G);
  const int k = (n - 3) / 6;
  const int a = n - 3;
  const int b = n - 2;
  const int c = n - 1;
  Graph T(n);
  add_edge(a, b, T);
  add_edge(b, c, T);
  add_edge(c, 4, T);
  for (int i = 0; i < k; ++i) {
    add_edge(6 * i + 0, 6 * i + 1, T);
    add_edge(6 * i + 1, 6 * i + 2, T);
    add_edge(6 * i + 1, 6 * i + 3, T);
    add_edge(6 * i + 1, 6 * i + 5, T);
    add_edge(6 * i + 4, 6 * i + 5, T);

    if (i < k - 1) {
      add_edge(6 * i + 1, 6 * (i + 1) + 4, T);
    }
  }
  return T;
}
