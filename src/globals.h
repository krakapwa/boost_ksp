#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <boost/graph/adjacency_list.hpp>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/shared_ptr.hpp>
#include <algorithm>    // std::find
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/transform_value_property_map.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/graph/copy.hpp>
#include <tuple>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>


namespace bp = boost::python;
namespace bn = boost::python::numpy;
namespace logging = boost::log;

using namespace boost;

namespace boost::trivial{
    enum logSeverityLevel
    {
        trace, //5
        debug, //4
        info, //3
        warning, //2
        error //1
    };
}

enum edge_myweight_t { edge_myweight };
enum edge_label_t { edge_label };
enum edge_id_t { edge_id };

struct GraphProperty{
    std::string name;
};

struct EdgeProperty{
    float weight;
    int label;
    int id;
    int id_vertex_in;
    int id_vertex_out;
};

struct VertexProperty{
    std::string name;
    int id;

    bool operator==(VertexProperty v){
        if(id == v.id)
            return true;
        else
            return false;
    }

    bool operator!=(VertexProperty v){
        if(id != v.id)
            return true;
        else
            return false;
    }
};

typedef adjacency_list<vecS, vecS, bidirectionalS,
    VertexProperty, EdgeProperty, GraphProperty> MyGraph;

typedef graph_traits<MyGraph>::vertex_iterator VertexIter;
typedef graph_traits<MyGraph>::edge_iterator EdgeIter;
typedef graph_traits<MyGraph>::edge_descriptor Edge;
typedef std::vector< graph_traits< MyGraph >::vertex_descriptor > VertexPath;
typedef std::vector< graph_traits< MyGraph >::edge_descriptor > EdgeSet;
typedef std::vector<EdgeSet> EdgeSets;
typedef std::vector<int> IdPath;
typedef graph_traits< MyGraph >::vertex_descriptor Vertex;
typedef std::tuple<EdgeSet, bool, std::vector<double>> ShortestPathRes;
#endif
