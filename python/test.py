
from graph_tool.all import *
from numpy.random import random
from random import choice, uniform
from numpy import sqrt, array, power, lexsort, transpose

def is_connected(g):
  _, histogram = label_components(g)
  return len(histogram) == 1

def connected_disk_graph(n, lo = 0.0, hi = 1.0, e = 1e-3):
  points = random((n, 2))
  points = points[lexsort(transpose(points)[::-1])]
  while True:
    r = (lo + hi) / 2
    g, pos = geometric_graph(points, r)
    if is_connected(g):
      if hi - lo < e:
        return g, pos
      else:
        hi = r
    else:
      if hi - lo < e / 1000:
        return None
      else:
        lo = r

class TreeVisitor(DFSVisitor):

    def __init__(self, tree_map):
        self.tree_map = tree_map

    def tree_edge(self, e):
        self.tree_map[e] = True

def dfs_tree(g):
  t = g.new_edge_property("bool")
  dfs_search(g, g.vertex(0), TreeVisitor(t))
  return t

def color_transform(a):
  return map(lambda x : 2 if x == 1 else 1 if x == 2 else 3, a)

def shape_transform(a):
  return map(lambda x : 2 if x == 1 else 0 if x == 2 else 6, a)

def weight_transform(a):
  return map(lambda x : 1.0 if x else 0.99, a)

def test(n, r = 0.0):
  g, p = connected_disk_graph(n, r)
  t = dfs_tree(g)

  g.set_edge_filter(t)

  fill = g.degree_property_map("total")
  fill.a = color_transform(fill.a)

  shape = g.degree_property_map("total")
  shape.a = shape_transform(shape.a)

  g.set_edge_filter(None)

  w = g.copy_property(t)
  w.a = weight_transform(t.a)

  pos = sfdp_layout(g, pos=p, eweight=w)

  pen = g.copy_property(t)
  pen.a = 2 * pen.a + 1.0

  control = g.new_edge_property("vector<double>")
  for e in g.edges():
    d = 0 if t[e] else choice([-1., 1.]) * uniform(0.3, 0.5) * sqrt(sum((pos[e.source()].a - pos[e.target()].a) ** 2))
    control[e] = [0.3, d, 0.7, d]

  graph_draw(g,
    pos=pos,
    vertex_fill_color=fill, vertex_shape=shape,
    edge_control_points=control,
    edge_color=t, edge_pen_width=pen)
