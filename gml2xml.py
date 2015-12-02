
import graph_tool.all as gt
from sys import argv
from re import findall

if __name__ == '__main__':
  for f in argv[1:]:
    g = gt.GraphView(gt.load_graph(f), directed=False, skip_properties=True)
    gt.remove_parallel_edges(g)
    gt.remove_self_loops(g)
    name = findall('[^/.]+', f)[-2].split('--')[0]
    g.save('output/{}.xml'.format(name))
    gt.graph_draw(g, output='output/{}.png'.format(name))

