
#include <boost/python.hpp>

#include "graph/bfs.hpp"
#include "graph/dfs.hpp"
#include "graph/rdfs.hpp"

BOOST_PYTHON_MODULE(graph)
{
  using namespace boost::python;
  def("bfs_tree", bfs_tree);
  def("dfs_tree", dfs_tree);
  def("rdfs_tree", rdfs_tree);
}

