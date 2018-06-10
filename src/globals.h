#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <boost/graph/adjacency_list.hpp>

namespace bp = boost::python;
namespace bn = boost::python::numpy;

using namespace boost;

enum edge_myweight_t { edge_myweight };
enum edge_label_t { edge_label };
enum edge_id_t { edge_id };

struct MyEdge{
    float weight;
    int label;
    int id;
};

typedef adjacency_list<vecS, vecS, bidirectionalS,
    no_property, MyEdge, no_property> MyGraph;

typedef graph_traits<MyGraph>::vertex_iterator VertexIter;
typedef graph_traits<MyGraph>::edge_iterator EdgeIter;
typedef graph_traits<MyGraph>::edge_descriptor Edge;
typedef std::vector< graph_traits< MyGraph >::vertex_descriptor > VertexPath;
typedef std::vector< graph_traits< MyGraph >::edge_descriptor > EdgePath;
typedef std::vector<int> IdPath;
typedef graph_traits< MyGraph >::vertex_descriptor Vertex;
#endif
