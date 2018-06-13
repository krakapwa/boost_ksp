#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/shared_ptr.hpp>
#include "pyksp.h"
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_internal_reference.hpp>

using namespace boost::python;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(add_edge_member_overloads,
                                       add_edge,
                                       3,
                                       7)

BOOST_PYTHON_MODULE(libksp)
{
  Py_Initialize();
  bn::initialize();

  class_<ksp, boost::shared_ptr<ksp> >("ksp",init<>())
    .def("create",&ksp::create )
    .staticmethod("create")
    .def("hello",&ksp::hello)
    .def("do_ksp",&ksp::do_ksp)
    .def("new_graph",&ksp::new_graph)
    .def("add_edge", &ksp::add_edge,
            add_edge_member_overloads(
              args("n0",
                   "n1",
                   "w",
                   "id",
                   "str_0",
                   "str_1",
                   "label"),
              "docstring add_edge")
    )
    .def("set_source",&ksp::set_source)
    .def("set_sink",&ksp::set_sink)
    .def("bellman_ford_shortest_paths",&ksp::bellman_ford_shortest_paths)
    ;
}
