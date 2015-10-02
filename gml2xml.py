
import graph_tool.all as gt
from sys import argv

if __name__ == '__main__':
  for f in argv[1:]:
    name = f.split('--')[0]
    g = gt.GraphView(gt.load_graph(f), directed=False, skip_properties=True)
    gt.remove_parallel_edges(g)
    gt.remove_self_loops(g)
    g.save('{}.xml'.format(name))

