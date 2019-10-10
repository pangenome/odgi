// odgi
#include "odgi.hpp"
//using namespace odgi;

// Pybind11
#include <pybind11/pybind11.h>
#include <pybind11/common.h>
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
//using namespace pybind11 as py;
namespace py = pybind11;

PYBIND11_MODULE(odgi, m)
{

    // Expose class Graph to Python.
    py::class_<odgi::graph_t>(m, "graph", "the odgi graph type")
        .def(py::init())
        .def("has_node",
              &odgi::graph_t::has_node,
             "Return true if the given node is in the graph.")
        .def("get_handle",
             &odgi::graph_t::get_handle,
             "Return the handle for the given node id.")
        .def("get_id",
             &odgi::graph_t::get_id,
             "Return the id of the given handle.",
             py::arg("handle"))
        .def("get_is_reverse",
             &odgi::graph_t::get_is_reverse,
             "Return true if the handle refers to the node reverse complement.",
             py::arg("handle"))
        .def("flip",
             &odgi::graph_t::flip,
             "Flip the handle to the opposite orientation.",
             py::arg("handle"))
        .def("get_length",
             &odgi::graph_t::get_length,
             "Return the length of the node referred to by the handle.",
             py::arg("handle"))
        .def("get_sequence",
             &odgi::graph_t::get_sequence,
             py::arg("handle"))
        .def("get_node_count",
             &odgi::graph_t::get_node_count,
             "Return the number of nodes in the graph.")
        .def("min_node_id",
             &odgi::graph_t::min_node_id,
             "Return the minimum node id in the graph.")
        .def("max_node_id",
             &odgi::graph_t::max_node_id,
             "Return the maximum node id in the graph.")
        // Definition of class_<odgi::graph_t> ends here.
    ;


}
