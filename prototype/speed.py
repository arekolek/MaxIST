

import graph_tool.all as gt


def read(f = int):
  return map(f, input().split())


class TreeVisitor(gt.DFSVisitor):

    def __init__(self, tree_map):
        self.tree_map = tree_map

    def tree_edge(self, e):
        self.tree_map[e] = True


def dfs_tree(g):
  t = g.new_edge_property("bool")
  gt.dfs_search(g, g.vertex(0), TreeVisitor(t))
  return gt.GraphView(g, efilt=t)


def rdfs_tree(g):
  return dfs_tree(g)


z, = read()
internal = 0
for i in range(z):
  n, m = read()
  g = gt.Graph()
  g.add_vertex(n)
  g.add_edge_list([tuple(read()) for j in range(m)])
  internal += (rdfs_tree(g).degree_property_map("total").a > 1).sum()
print(internal / z)

