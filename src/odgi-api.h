// odgi-api.h ODGI C-API for FFIs
//
// Copyright (c) 2022 Erik Garrison and Pjotr Prins
//
// ODGI is published under the MIT license

#pragma once

#include <cstddef>

#include "odgi.hpp"

extern "C" {

  using namespace odgi;

  // Introduce opaque types to support type checking of pointers to C++ classes
  typedef struct opaque_graph {} *ograph_t;

  // These functions are exposed as a simple C API
  const char *odgi_version();
  const ograph_t odgi_load_graph(const char *filen);
  void odgi_free_graph(const ograph_t graph);
  const size_t odgi_get_node_count(ograph_t graph);
  const size_t odgi_max_node_id(const ograph_t graph);
  const size_t odgi_min_node_id(const ograph_t graph);
  const size_t odgi_get_path_count(ograph_t graph);
  void odgi_for_each_path_handle(const ograph_t graph, void (*next) (const path_handle_t path));
  const bool odgi_for_each_handle(const ograph_t graph,
                                  bool (*next) (const handle_t handle));
  const bool odgi_follow_edges(const ograph_t graph,
                               const handle_t handle,
                               bool go_left,
                               bool (*next) (const handle_t handle));
  const handle_t odgi_edge_first_handle(const ograph_t graph, const edge_t &edge_handle);
  const handle_t odgi_edge_second_handle(const ograph_t graph, const edge_t &edge_handle);
  const bool odgi_has_node(const ograph_t graph, nid_t node_id);
  const char *odgi_get_sequence(const ograph_t graph, const handle_t handle);
  const nid_t odgi_get_id(const ograph_t graph, const handle_t handle);
  const bool odgi_get_is_reverse(const ograph_t graph, const handle_t handle);
  const size_t odgi_get_length(const ograph_t graph, const handle_t handle);
  // Path handling
  const char *odgi_get_path_name(const ograph_t graph, const path_handle_t path);
  const bool odgi_has_path(const ograph_t graph, const char *path_name);
  const bool odgi_path_is_empty(const ograph_t graph, const path_handle_t path);
  const path_handle_t odgi_get_path_handle(const ograph_t graph, const char *path_name);
  // Steps
  const size_t odgi_get_step_count(const ograph_t graph, const handle_t handle);
  const handle_t odgi_step_get_handle(const ograph_t graph, step_handle_t step);
  const path_handle_t odgi_step_get_path(const ograph_t graph, step_handle_t step);
  const step_handle_t odgi_step_path_begin(const ograph_t graph, path_handle_t path);
  const step_handle_t odgi_step_path_end(const ograph_t graph, path_handle_t path);
  const step_handle_t odgi_step_path_back(const ograph_t graph, path_handle_t path);
  const int64_t odgi_step_path_id(const ograph_t graph, step_handle_t &step_handle);
  bool odgi_step_is_reverse(const ograph_t graph, step_handle_t &step_handle);
  const int64_t odgi_step_prev_id(const ograph_t graph, step_handle_t &step_handle);
  const int64_t odgi_step_prev_rank(const ograph_t graph, step_handle_t &step_handle);
  const int64_t odgi_step_next_id(const ograph_t graph, step_handle_t &step_handle);
  const int64_t odgi_step_next_rank(const ograph_t graph,step_handle_t &step_handle);
  const bool odgi_step_eq(const ograph_t graph,
                          step_handle_t &step_handle1,
                          step_handle_t &step_handle2);
  const step_handle_t odgi_path_front_end(const ograph_t graph,
                                     const path_handle_t path);
  const step_handle_t odgi_get_next_step(const ograph_t graph,
                                         const step_handle_t step);
  const step_handle_t odgi_get_previous_step(const ograph_t graph,
                                         const step_handle_t step);
  const bool odgi_has_edge(const ograph_t graph, const handle_t left, const handle_t right);
  const bool odgi_is_path_front_end(const ograph_t graph, const step_handle_t step);
  const bool odgi_is_path_end(const ograph_t graph, const step_handle_t step);
  const bool odgi_has_next_step(const ograph_t graph, const step_handle_t step);
  const bool odgi_has_previous_step(const ograph_t graph, const step_handle_t step);
  const path_handle_t odgi_get_path_handle_of_step(const ograph_t graph,
                                                   const step_handle_t step);
  void odgi_for_each_step_in_path(const ograph_t graph,
                                  const path_handle_t path,
                                  void (*next) (const step_handle_t step));
  const bool odgi_for_each_step_on_handle(const ograph_t graph,
                                          const handle_t handle,
                                          bool (*next) (const step_handle_t step));

} // extern "C"
