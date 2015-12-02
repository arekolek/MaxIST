
import graph_tool.all as gt

from numpy import log, linalg, transpose, triu, sqrt
from numpy.random import random, permutation


def graph_union(g, h):
  g = gt.Graph(g, prune=True)
  in_h_notin_g = (gt.adjacency(h) - gt.adjacency(g)).toarray() > 0
  edges = transpose(triu(in_h_notin_g, 1).nonzero())
  g.add_edge_list(edges)
  return g


def random_graph(n, p):
	g = gt.Graph(directed=False)
	g.add_vertex(n)
	g.add_edge_list(transpose(triu(random((n, n)) < p, 1).nonzero()))
	return g


def connected_random(n, p):
  g = gt.complete_graph(n)
  weight = g.new_edge_property("double", random(g.num_edges()))
  t = gt.GraphView(g, efilt=gt.min_spanning_tree(g, weights=weight))
  h = gt.GraphView(g, efilt=g.new_edge_property("bool", weight.a < p))
  return graph_union(h, t)


def comparison(n, pp, pt):
  g = gt.complete_graph(n)
  weight = g.new_edge_property("double", random(g.num_edges()))
  
  t = gt.min_spanning_tree(g, weights=weight)
  
  p = g.new_edge_property('bool')
  q = permutation(n)
  for i in range(n-1):
    p[g.edge(q[i], q[i+1])] = True
  
  hp = g.new_edge_property("bool", weight.a < pp)
  ht = g.new_edge_property("bool", weight.a < pt)
  
  draw2(g, hp, p, 'gnp+path-{}.png'.format(pp))
  draw2(g, ht, p, 'gnp+path-{}.png'.format(pt))
  draw2(g, hp, t, 'gnp+mst-{}.png'.format(pp))
  draw2(g, ht, t, 'gnp+mst-{}.png'.format(pt))


def comparison2(n, rp, rt):
  g = gt.complete_graph(n)
  pts = random((n, 2))
  
  p = g.new_edge_property('bool')
  q = permutation(n)
  for i in range(n-1):
    p[g.edge(q[i], q[i+1])] = True
  
  d, pos = gt.triangulation(pts, type="delaunay")
  t = create_filter(g, gt.GraphView(d, efilt=min_spanning_tree(d, pos)))
  
  hp = create_filter(g, gt.geometric_graph(pts, rp)[0])
  ht = create_filter(g, gt.geometric_graph(pts, rt)[0])
  
  draw2(g, hp, p, 'rgg+path-{}.png'.format(rp))
  draw2(g, ht, p, 'rgg+path-{}.png'.format(rt))
  draw2(g, ht, t, 'rgg+mst-{}.png'.format(rt))

    
def create_filter(g, h):
  f = g.new_edge_property('bool')
  for e in h.edges():
    f[g.edge(h.vertex_index[e.source()], h.vertex_index[e.target()])] = True
  return f
  

def layout(g, h, t, pos=None):
  g = gt.GraphView(g, efilt=lambda e: h[e] or t[e])
  return gt.sfdp_layout(g, pos=pos)

def draw2(g, h, t, output, pos=None):
  g = gt.GraphView(g, efilt=lambda e: h[e] or t[e])
  
  print('{} {}'.format(output, 2*g.num_edges()/g.num_vertices()))

  width = g.new_edge_property('double')
  color = g.new_edge_property('vector<double>')
  
  for e in g.edges():
    width[e] = 1.6 if h[e] else 0.8
    color[e] = (0.1,0.9,0.1,1) if h[e] and t[e] else (0.9,0.1,0.1,1) if t[e] else (0,0,0,1)
  
  gt.graph_draw(g, pos=pos,
    vertex_size=8, vertex_color='black', vertex_fill_color='black',
    edge_pen_width=1.6, edge_color=color,
    output_size=(500, 500), output=output, bg_color=[1,1,1,1])

def random_geometric(n, r):
	return gt.geometric_graph(random((n, 2)), r)


def draw(g, output='graph.png', pos=None):
	gt.graph_draw(g, pos=pos, vertex_size=8, vertex_color='black', vertex_fill_color='black', edge_pen_width=0.8, edge_color='black',
			   output_size=(500, 500), output=output, bg_color=[1,1,1,1])


def threshold(n, e, c = 1):
	p = c*pow(n, e)
	draw(random_graph(n, p), 'gnp_100_{:0>6}.png'.format(int(1000000*p)))
	print(p)


def gnp():
	n = 100
	p = [#(-2-.1, 0.5), # empty graphs
		 #(-2+.1, 3), # non-empty graphs
		 #(-1.5+.1, 1.7), # edges have a vertex in common
		 (-1-1/3+.1, 1.4), # trees on 3+1 vertices
		 (-1+.1, 0.85), # triangles and cycles
		 (-1, 1.1*log(n)), # connected
		 #(-2/3, 1.6), # K4 subgraphs
		 #(-2/5, 1), # K6 subgraphs
		 #(-0.5, sqrt(log(n))),  # vertices have a common neighbor
		 ]

	for e, c in p:
		threshold(n, e, c)


def add_random_path(h):
  g = gt.Graph(directed=False)
  n = h.num_vertices()
  q = permutation(n)
  g.add_edge_list([(q[i], q[i+1]) for i in range(n-1)])
  return graph_union(h, g)


def gnp(n, p, path=False, isolated=True):
	g = random_graph(n, p)
	if path:
		add_random_path(g)
	if not isolated:
		g.set_vertex_filter(g.new_vertex_property('bool', g.degree_property_map('total').a > 0))
	draw(g, 'gnp_{}_{:0>3}.png'.format(n, int(1000*p)))


def rgg(n, p, path=False, isolated=True):
	g = random_geometric(n, p)[0]
	if path:
		add_random_path(g)
	if not isolated:
		g.set_vertex_filter(g.new_vertex_property('bool', g.degree_property_map('total').a > 0))
	draw(g, 'rgg_{}_{:0>3}.png'.format(n, int(1000*p)))


def min_spanning_tree(g, pos):
  weight = g.new_edge_property("double")
  for e in g.edges():
    weight[e] = linalg.norm(pos[e.target()].a - pos[e.source()].a)
  return gt.min_spanning_tree(g, weights=weight)


def connected_geometric(n, r):
  pts = random((n, 2))
  g, pos = gt.triangulation(pts, type="delaunay")
  t = gt.GraphView(g, efilt=min_spanning_tree(g, pos))
  h, pos = gt.geometric_graph(pts, r)
  h = graph_union(h, t)
  return h, pos


if __name__ == '__main__':
	gnp()
	for x in range(5, 16, 5):
		rgg(x/100)
	
