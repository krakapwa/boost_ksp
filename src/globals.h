#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <boost/graph/adjacency_list.hpp>

namespace bp = boost::python;
namespace bn = boost::python::numpy;

using namespace boost;

enum edge_myweight_t { edge_myweight };
enum edge_label_t { edge_label };

namespace boost {
  BOOST_INSTALL_PROPERTY(edge, myweight);
  BOOST_INSTALL_PROPERTY(edge, label);
}

typedef property<edge_label_t, int> Label;
typedef property<edge_myweight_t, float, Label> Weight;
typedef adjacency_list<vecS, vecS, bidirectionalS,
    no_property, Weight> MyGraph;
typedef graph_traits<MyGraph>::vertex_iterator vertex_iter;
typedef std::vector< graph_traits< MyGraph >::vertex_descriptor > Path;
typedef graph_traits< MyGraph >::vertex_descriptor Vertex;
typedef property_map<MyGraph, vertex_index_t>::type IndexMap;
#endif
