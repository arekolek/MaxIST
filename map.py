from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

from xml.etree import ElementTree as ET
from itertools import chain
import graph_tool.all as gt


def vertices(f):
  root = ET.parse(f)
  ns = {'s': 'http://sndlib.zib.de/network'}
  for e in root.findall('*/*/s:node', ns):
    yield (float(e.find('*/s:x', ns).text), float(e.find('*/s:y', ns).text)) 


def edges(f):
  root = ET.parse(f)
  ns = {'s': 'http://sndlib.zib.de/network'}
  pos = {e.get('id'):(float(e.find('*/s:x', ns).text), float(e.find('*/s:y', ns).text)) for e in root.findall('*/*/s:node', ns)}
  for e in root.findall('*/*/s:link', ns):
    yield chain(pos[e.find('s:source', ns).text], pos[e.find('s:target', ns).text])


if __name__ == '__main__':
  from sys import argv
  from re import findall
  
  for f in argv[1:]:
    vs = np.array(list(vertices(f)))

    xmin, ymin = vs.min(axis=0)
    xmax, ymax = vs.max(axis=0)
    
    #x, y = vs.mean(axis=0)
    x, y = (xmin+xmax)/2, (ymin+ymax)/2
    
    m = Basemap(projection='stere', lon_0=x, lat_0=y, width=1000, height=1000)
    xlo, ylo = m(xmin, ymin)
    xhi, yhi = m(xmax, ymax)
    span = max(xhi-xlo, yhi-ylo) * 1.15
    
    #xmin, xmax = xmin-(xmax-xmin)/10, xmax+(xmax-xmin)/10
    #ymin, ymax = ymin-(ymax-ymin)/5, ymax+(ymax-ymin)/5

    # create new figure, axes instances.
    fig = plt.figure(frameon=False)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    # setup mercator map projection.
    m = Basemap(#llcrnrlon=xmin,llcrnrlat=ymin,urcrnrlon=xmax,urcrnrlat=ymax,\
                #rsphere=(6378137.00,6356752.3142),\
                resolution='l', projection='stere',\
                lon_0=x, lat_0=y, width=span, height=span)

    #m.drawcountries(linestyle='dotted')
    m.fillcontinents(color='#dddddd')
    
    for e in edges(f):
      m.drawgreatcircle(*e,linewidth=0.2,color='black')
    
    m.scatter(vs[:,0], vs[:,1], latlon=True, s=4, c='black', alpha=1, zorder=10)
    
    name = findall('[^/.]+', f)[-2]
    fig.set_size_inches(1.042, 1.042)
    fig.savefig('output/{}.pdf'.format(name), dpi=100)
    
