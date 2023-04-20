// odgi ffi calls C-API functions
//
// Copyright © 2022 Pjotr Prins

// Pybind11
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/iostream.h>

#include "odgi.hpp"
#include "odgi-api.h"

PYBIND11_MAKE_OPAQUE(ograph_t);

namespace py = pybind11;

PYBIND11_MODULE(odgi_ffi, m)
{
    py::class_<ograph_t>(m, "opaque_graph pointer for FFI");
    // py::class_<path_handle_t>(m, "opaque path_handle_t for FFI");
    // py::class_<handle_t>(m, "opaque handle_t for FFI");
    py::class_<step_handle_t>(m, "step_handle_t for FFI");
    m.def("odgi_version", &odgi_version, "Get the odgi library build version");
    m.def("odgi_long_long_size", &odgi_long_long_size);
    m.def("odgi_handle_i_size", &odgi_handle_i_size);
    m.def("odgi_step_handle_i_size", &odgi_step_handle_i_size);
    m.def("odgi_test_uint128", &odgi_test_uint128);
    m.def("odgi_load_graph", &odgi_load_graph);
    m.def("odgi_free_graph", &odgi_free_graph);
    m.def("odgi_get_node_count", &odgi_get_node_count);
    m.def("odgi_max_node_id", &odgi_max_node_id);
    m.def("odgi_min_node_id", &odgi_min_node_id);
    m.def("odgi_get_path_count", &odgi_get_path_count);
    m.def("odgi_for_each_path_handle", &odgi_for_each_path_handle);
    m.def("odgi_for_each_handle", &odgi_for_each_handle);
    m.def("odgi_follow_edges", &odgi_follow_edges);
    m.def("odgi_edge_first_handle", &odgi_edge_first_handle);
    m.def("odgi_edge_second_handle", &odgi_edge_second_handle);
    m.def("odgi_has_node", &odgi_has_node);
    m.def("odgi_get_sequence", &odgi_get_sequence);
    m.def("odgi_get_id", &odgi_get_id);
    m.def("odgi_get_is_reverse", &odgi_get_is_reverse);
    m.def("odgi_get_length", &odgi_get_length);
    m.def("odgi_get_path_name", &odgi_get_path_name);
    m.def("odgi_has_path", &odgi_has_path);
    m.def("odgi_path_is_empty", &odgi_path_is_empty);
    m.def("odgi_get_path_handle", &odgi_get_path_handle);
    m.def("odgi_get_step_count", &odgi_get_step_count);
    m.def("odgi_get_handle_of_step", &odgi_get_handle_of_step);
    m.def("odgi_get_path", &odgi_get_path);
    m.def("odgi_path_begin", &odgi_path_begin);
    m.def("odgi_path_end", &odgi_path_end);
    m.def("odgi_path_back", &odgi_path_back);
    m.def("odgi_step_path_id", &odgi_step_path_id);
    m.def("odgi_step_is_reverse", &odgi_step_is_reverse);
    m.def("odgi_step_prev_id", &odgi_step_prev_id);
    m.def("odgi_step_prev_rank", &odgi_step_prev_rank);
    m.def("odgi_step_next_id", &odgi_step_next_id);
    m.def("odgi_step_next_rank", &odgi_step_next_rank);
    m.def("odgi_step_eq", &odgi_step_eq);
    m.def("odgi_path_front_end", &odgi_path_front_end);
    m.def("odgi_get_next_step", &odgi_get_next_step);
    m.def("odgi_get_previous_step", &odgi_get_previous_step);
    m.def("odgi_has_edge", &odgi_has_edge);
    m.def("odgi_is_path_front_end", &odgi_is_path_front_end);
    m.def("odgi_is_path_end", &odgi_is_path_end);
    m.def("odgi_has_next_step", &odgi_has_next_step);
    m.def("odgi_has_previous_step", &odgi_has_previous_step);
    m.def("odgi_get_path_handle_of_step", &odgi_get_path_handle_of_step);
    m.def("odgi_for_each_step_in_path", &odgi_for_each_step_in_path);
    m.def("odgi_for_each_step_on_handle", &odgi_for_each_step_on_handle);
}
