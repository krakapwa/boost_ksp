#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_internal_reference.hpp>
#include "ksp.h"

using namespace boost::python;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(add_edge_member_overloads,
                                       add_edge,
                                       3,
                                       7)

BOOST_PYTHON_MODULE(libksp)
{
  Py_Initialize();
  namespace bp = boost::python;

  class_<Ksp>("ksp")
    .def("run",&Ksp::run)
    .def("add_edge", &Ksp::add_edge,
            add_edge_member_overloads(
              args("n0",
                   "n1",
                   "w",
                   "id",
                   "str_0",
                   "str_1",
                   "label"),
              "docstring add_edge"))
    .def("config", &Ksp::config,
         (bp::arg("source_vertex_id"),
          bp::arg("sink_vertex_id"),
          bp::arg("l_max")=-1,
          bp::arg("source_vertex_name")="source",
          bp::arg("sink_vertex_name")="sink",
          bp::arg("loglevel")="info",
          bp::arg("min_cost")=false))
    .def("set_loglevel",&Ksp::set_loglevel)
    .def("set_source",&Ksp::set_source)
    .def("set_sink",&Ksp::set_sink)
    .def("remove_edge",&Ksp::remove_edge)
    .def("set_label_all_edges",&Ksp::set_label_all_edges)
    .def("num_edges",&Ksp::num_edges)
    .def("num_vertices",&Ksp::num_vertices)
    ;
}
