#include <lost.hpp>
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "LOST tests"
#include <boost/test/unit_test.hpp>

boost::adjacency_list<boost::hash_setS, boost::vecS, boost::undirectedS> typedef graph;

BOOST_AUTO_TEST_CASE( lost_test )
{
  // 0-1-2-3-4-5-6
  //     |   |
  //     7   9
  //     |   |
  //     8   10
  graph T;
  add_edge(0, 1, T);
  add_edge(1, 2, T);
  add_edge(2, 3, T);
  add_edge(3, 4, T);
  add_edge(4, 5, T);
  add_edge(5, 6, T);
  add_edge(2, 7, T);
  add_edge(7, 8, T);
  add_edge(4, 9, T);
  add_edge(9, 10, T);
  leaf_info<graph> info(T);
  BOOST_CHECK( !info.is_path() );
  for(int i = 0; i < 3; ++i)
    BOOST_CHECK( info.on_branch(0, i) );
  for(int i = 3; i < 11; ++i)
    BOOST_CHECK( !info.on_branch(0, i) );
}
