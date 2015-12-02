
import graph_tool.all as gt
import numpy as np


def scale(rawpoints, low=0.0, high=100.0):
    mins = np.min(rawpoints, axis=0)
    maxs = np.max(rawpoints, axis=0)
    rng = maxs - mins
    return high - (((high - low) * (maxs - rawpoints)) / rng)


def sndlib(f):
  from xml.etree import ElementTree as ET

  root = ET.parse(f)
  ns = {'s': 'http://sndlib.zib.de/network'}
  
  g = gt.Graph(directed=False)
  g.add_vertex(len(root.findall('*/*/s:node', ns)))
  
  index = {e.get('id'):i for i,e in enumerate(root.findall('*/*/s:node', ns))}
  for e in root.findall('*/*/s:link', ns):
    g.add_edge(index[e.find('s:source', ns).text], index[e.find('s:target', ns).text])
  pos = g.new_vertex_property('vector<float>', scale(np.array([(float(e.find('*/s:x', ns).text), -float(e.find('*/s:y', ns).text)) for e in root.findall('*/*/s:node', ns)])))
  
  gt.remove_parallel_edges(g)
  gt.remove_self_loops(g)
  
  return g, pos


if __name__ == '__main__':
  from sys import argv
  from re import findall
  from statistics import mean
  
  for f in argv[1:]:
    g, pos = sndlib(f)
    
    name = findall('[^/.]+', f)[-2]
    
    deg = g.degree_property_map('total').a
    print(name, g.num_vertices(), g.num_edges(), '{:.2f}'.format(mean(deg)), min(deg), max(deg), sep=' & ', end=' \\\\\n')
    
    g.save('output/{}.xml'.format(name))
    gt.graph_draw(g, pos, vertex_size=2, vertex_color='black', vertex_fill_color='black', edge_pen_width=0.2, edge_color='black', output='output/{}.pdf'.format(name), output_size=(75,75))

