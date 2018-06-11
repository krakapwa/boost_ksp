#ifndef PYKSP_H
#define PYKSP_H

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include "boost_ksp.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/transform_value_property_map.hpp>
#include "globals.h"
#include <tuple>

namespace bp = boost::python;
namespace bn = boost::python::numpy;
using namespace boost;


struct ksp {

    MyGraph * G;

    // gets the weight property
    //property_map<MyGraph, edge_myweight_t>::type weight_pmap =
    //    get(edge_myweight_t(), *G);

    // gets the label property
    //property_map<MyGraph, edge_label_t>::type label_pmap =
    //    get(edge_label_t(), *G);

    // gets the id property
    //property_map<MyGraph, edge_id_t>::type id_pmap =
    //    get(edge_id_t(), *G);

    int n_vertices;
    Vertex source_vertex;
    Vertex sink_vertex;


    static shared_ptr<ksp> create();
    void new_graph(int n_vertices);

    bool do_ksp();
    void print_edges();
    Vertex my_add_vertex(int id, std::string str);
    void add_edge(int n0, int n1, double w, int id,
                  std::string str_0, std::string str_1);
    void invert_edge(Edge e, bool inv_algebraic_sign);
    void invert_edges(EdgePath edge_path, bool inv_algebraic_sign);
    void set_source(int id, std::string str);
    void set_sink(int id, std::string str);
    ShortestPathRes bellman_ford_shortest_paths();
    ShortestPathRes dijkstra_shortest_paths();

    void cost_transform(const std::vector<double> & distance);

    std::string hello() { return "Just nod if you can hear me!"; }

    void run(int n0,
                    int n1,
                    double w);
};

#endif
