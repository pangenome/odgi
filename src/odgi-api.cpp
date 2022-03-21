// odgi-api exports the C-API for FFIs
//
// Copyright (c) 2022 Erik Garrison and Pjotr Prins
//
// ODGI is published under the MIT license
//
// IMPORTANT: this interface is WIP. Function names reflect odgi
// naming and we suggest to keep that as a low-level interface. A more
// natural OOP representation may be built on top of this low-level
// API.
//
// The current edition is in memory graph read-only. A graph write
// edition is future work.

#include "odgi-api.h"
#include "version.hpp"

extern "C" {

  using namespace odgi;

  const char *odgi_version() {
    return Version::get_version().c_str();
  }

  const ograph_t odgi_load_graph(const char *filen) {
    graph_t *graph = new graph_t();
    std::ifstream in(filen);
    graph->deserialize(in);
    return (ograph_t)graph;
  }

  // Calling this function is probably not safe
  void odgi_free_graph(ograph_t graph) {
    graph_t *g = (graph_t *)graph;
    delete g;
  }

  const size_t odgi_get_node_count(ograph_t graph) {
    return ((graph_t *)graph)->get_node_count();
  }

  const size_t odgi_max_node_id(const ograph_t graph) {
    return ((graph_t *)graph)->max_node_id();
  }

  const size_t odgi_min_node_id(const ograph_t graph) {
    return ((graph_t *)graph)->min_node_id();
  }

  const size_t odgi_get_path_count(ograph_t graph) {
    return ((graph_t *)graph)->get_path_count();
  }

  void odgi_for_each_path_handle(const ograph_t graph, void (*next) (const path_handle_t path)) {
    ((graph_t *)graph)->for_each_path_handle([&](const path_handle_t& path) {
      next(path);
    });
  }

  const bool odgi_for_each_handle(const ograph_t graph,
                                  // const std::function<bool(const handle_t&)>& iteratee)
                                  bool (*next) (const handle_t handle))
  {
    return ((graph_t *)graph)->for_each_handle([&](const handle_t h)
    {
      next(h);
      return true;
    },
      false); // note that parallel is false
  }

  /// Loop over all the handles to next/previous (right/left) nodes. Passes
  /// them to a callback which returns false to stop iterating and true to
  /// continue. Returns true if we finished and false if we stopped early.
  const bool odgi_follow_edges(const ograph_t graph,
                               const handle_t handle,
                               bool go_left,
                               bool (*next) (const handle_t handle))
  {
    return ((graph_t *)graph)->follow_edges(handle,go_left,
                               [&](const handle_t handle)
                               {
                                 next(handle);
                                 return true;
                               }
                               );
  };

  const handle_t odgi_edge_first_handle(const ograph_t graph, const edge_t &edge_handle)
  {
      return (&edge_handle)->first;
  }

  const handle_t odgi_edge_second_handle(const ograph_t graph, const edge_t &edge_handle)
  {
      return (&edge_handle)->second;
  }

  const bool odgi_has_node(const ograph_t graph, nid_t node_id) {
    return ((graph_t *)graph)->has_node(node_id);
  }

  const char *odgi_get_sequence(const ograph_t graph, const handle_t handle) {
    return ((graph_t *)graph)->get_sequence(handle).c_str();
  }

  const nid_t odgi_get_id(const ograph_t graph, const handle_t handle) {
    return ((graph_t *)graph)->get_id(handle);
  }

  const bool odgi_get_is_reverse(const ograph_t graph, const handle_t handle) {
    return ((graph_t *)graph)->get_is_reverse(handle);
  }

  const size_t odgi_get_length(const ograph_t graph, const handle_t handle) {
    return ((graph_t *)graph)->get_length(handle);
  }

  // Path handling
  const char *odgi_get_path_name(const ograph_t graph, const path_handle_t path) {
    return ((graph_t *)graph)->get_path_name(path).c_str();
  }

  const bool odgi_has_path(const ograph_t graph, const char *path_name) {
    return ((graph_t *)graph)->has_path(path_name);
  }

  const bool odgi_path_is_empty(const ograph_t graph, const path_handle_t path) {
    return ((graph_t *)graph)->is_empty(path);
  }

  const path_handle_t odgi_get_path_handle(const ograph_t graph, const char *path_name) {
    return ((graph_t *)graph)->get_path_handle(path_name);
  }

  // Steps
  const size_t odgi_get_step_count(const ograph_t graph, const handle_t handle)
  {
    return ((graph_t *)graph)->get_step_count(handle);
  }

  const handle_t odgi_step_get_handle(const ograph_t graph, step_handle_t step)
  {
    return ((graph_t *)graph)->get_handle_of_step(step);
  }

  const path_handle_t odgi_step_get_path(const ograph_t graph, step_handle_t step)
  {
    return ((graph_t *)graph)->get_path(step);
  }

  const step_handle_t odgi_step_path_begin(const ograph_t graph, path_handle_t path)
  {
    return ((graph_t *)graph)->path_begin(path);
  }

  const step_handle_t odgi_step_path_end(const ograph_t graph, path_handle_t path)
  {
    return ((graph_t *)graph)->path_end(path);
  }

  const step_handle_t odgi_step_path_back(const ograph_t graph, path_handle_t path)
  {
    return ((graph_t *)graph)->path_back(path);
  }

  const int64_t odgi_step_path_id(const ograph_t graph, step_handle_t &step_handle) {
    return (reinterpret_cast<const int64_t*>(&step_handle)[0]) >> 1;
  }

  bool odgi_step_is_reverse(const ograph_t graph, step_handle_t &step_handle) {
    return (reinterpret_cast<const int64_t*>(&step_handle)[0]) & 1;
  }

  const int64_t odgi_step_prev_id(const ograph_t graph, step_handle_t &step_handle)
  {
      return (reinterpret_cast<const int64_t*>(&step_handle)[1]);
  }

  const int64_t odgi_step_prev_rank(const ograph_t graph, step_handle_t &step_handle)
  {
      return (reinterpret_cast<const int64_t*>(&step_handle)[2]);
  }

  const int64_t odgi_step_next_id(const ograph_t graph, step_handle_t &step_handle)
  {
      return (reinterpret_cast<const int64_t*>(&step_handle)[3]);
  }

  const int64_t odgi_step_next_rank(const ograph_t graph,step_handle_t &step_handle)
  {
      return (reinterpret_cast<const int64_t*>(&step_handle)[4]);
  }

  const bool odgi_step_eq(const ograph_t graph,
                          step_handle_t &step_handle1,
                          step_handle_t &step_handle2)
  {
    // this code is hardly optimal. FIXME.
    for (uint i=0;i<8;i++) {
      if ((reinterpret_cast<const int64_t*>(&step_handle1)[i]) != (reinterpret_cast<const int64_t*>(&step_handle2)[i])){
        return false;
      }
    }
    return true;
  }

  const step_handle_t odgi_path_front_end(const ograph_t graph,
                                     const path_handle_t path)
  {
    return ((graph_t *)graph)->path_front_end(path);
  }

  const step_handle_t odgi_get_next_step(const ograph_t graph,
                                         const step_handle_t step)
  {
    return ((graph_t *)graph)->get_next_step(step);
  }

  const step_handle_t odgi_get_previous_step(const ograph_t graph,
                                         const step_handle_t step)
  {
    return ((graph_t *)graph)->get_previous_step(step);
  }

  const bool odgi_has_edge(const ograph_t graph, const handle_t left, const handle_t right)
  {
    return ((graph_t *)graph)->has_edge(left,right);
  }

  const bool odgi_is_path_front_end(const ograph_t graph, const step_handle_t step)
  {
    return ((graph_t *)graph)->is_path_front_end(step);
  }
  const bool odgi_is_path_end(const ograph_t graph, const step_handle_t step)
  {
    return ((graph_t *)graph)->is_path_end(step);
  }
  const bool odgi_has_next_step(const ograph_t graph, const step_handle_t step)
  {
    return ((graph_t *)graph)->has_next_step(step);
  }
  const bool odgi_has_previous_step(const ograph_t graph, const step_handle_t step)
  {
    return ((graph_t *)graph)->has_previous_step(step);
  }

  const path_handle_t odgi_get_path_handle_of_step(const ograph_t graph,
                                                   const step_handle_t step)
  {
    return ((graph_t *)graph)->get_path_handle_of_step(step);
  }

  void odgi_for_each_step_in_path(const ograph_t graph,
                                  const path_handle_t path,
                                  void (*next) (const step_handle_t step))
  {
    ((graph_t *)graph)->for_each_step_in_path(path,
                                 [&](const step_handle_t step) {
                                   next(step);
                                 });
  };


  const bool odgi_for_each_step_on_handle(const ograph_t graph,
                                          const handle_t handle,
                                          bool (*next) (const step_handle_t step))
  {
    return ((graph_t *)graph)->for_each_step_on_handle(handle,
                                          [&](const step_handle_t step) {
                                            next(step);
                                            return true; });
  }

} // extern "C"
