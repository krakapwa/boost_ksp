#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/shared_ptr.hpp>
#include "pyksp.h"

namespace bp = boost::python;
namespace bn = boost::python::numpy;
using namespace boost;

BOOST_PYTHON_MODULE(libksp)
{
  using namespace boost::python;
  Py_Initialize();
  bn::initialize();

  class_<ksp, boost::shared_ptr<ksp> >("ksp",init<>())
    .def("create",&ksp::create )
    .staticmethod("create")
    .def("hello",&ksp::hello)
    .def("new_graph",&ksp::new_graph)
    .def("add_edge",&ksp::add_edge)
    .def("set_source_id",&ksp::set_source_id)
    .def("set_sink_id",&ksp::set_sink_id)
    .def("bellman_ford_shortest_paths",&ksp::bellman_ford_shortest_paths)
    ;
}
