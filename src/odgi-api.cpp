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

  const graph_t* odgi_load_graph(const char *filen) {
    graph_t *graph = new graph_t();
    std::ifstream in(filen);
    graph->deserialize(in);
    return graph;
  }

  void odgi_free_graph(graph_t* graph) {
    delete graph;
  }

  const size_t odgi_get_node_count(const graph_t* graph) {
    return graph->get_node_count();
  }

  const size_t odgi_max_node_id(const graph_t* graph) {
    return graph->max_node_id();
  }

  const size_t odgi_min_node_id(const graph_t* graph) {
    return graph->min_node_id();
  }

  const size_t odgi_get_path_count(const graph_t* graph) {
    return graph->get_path_count();
  }

  void odgi_for_each_path_handle(const graph_t *graph, void (*next) (const path_handle_t path)) {
    graph->for_each_path_handle([&](const path_handle_t& path) {
      next(path);
    });
  }

  const bool odgi_for_each_handle(const graph_t *graph,
                                  // const std::function<bool(const handle_t&)>& iteratee)
                                  bool (*next) (const handle_t handle))
  {
    return graph->for_each_handle([&](const handle_t h)
    {
      next(h);
      return true;
    },
      false); // note that parallel is false
  }

  /// Loop over all the handles to next/previous (right/left) nodes. Passes
  /// them to a callback which returns false to stop iterating and true to
  /// continue. Returns true if we finished and false if we stopped early.
  const bool odgi_follow_edges(const graph_t *graph,
                               const handle_t handle,
                               bool go_left,
                               bool (*next) (const handle_t handle))
  {
    return graph->follow_edges(handle,go_left,
                               [&](const handle_t handle)
                               {
                                 next(handle);
                                 return true;
                               }
                               );
  };

  const handle_t odgi_edge_first_handle(const graph_t *graph, const edge_t &edge_handle)
  {
      return (&edge_handle)->first;
  }

  const handle_t odgi_edge_second_handle(const graph_t *graph, const edge_t &edge_handle)
  {
      return (&edge_handle)->second;
  }

  const bool odgi_has_node(const graph_t *graph, nid_t node_id) {
    return graph->has_node(node_id);
  }

  const char *odgi_get_sequence(const graph_t *graph, const handle_t handle) {
    return graph->get_sequence(handle).c_str();
  }

  const nid_t odgi_get_id(const graph_t *graph, const handle_t handle) {
    return graph->get_id(handle);
  }

  const bool odgi_get_is_reverse(const graph_t *graph, const handle_t handle) {
    return graph->get_is_reverse(handle);
  }

  const size_t odgi_get_length(const graph_t *graph, const handle_t handle) {
    return graph->get_length(handle);
  }

  // Path handling
  const char *odgi_get_path_name(const graph_t *graph, const path_handle_t path) {
    return graph->get_path_name(path).c_str();
  }

  const bool odgi_has_path(const graph_t *graph, const char *path_name) {
    return graph->has_path(path_name);
  }

  const bool odgi_path_is_empty(const graph_t *graph, const path_handle_t path) {
    return graph->is_empty(path);
  }

  const path_handle_t odgi_get_path_handle(const graph_t *graph, const char *path_name) {
    return graph->get_path_handle(path_name);
  }

  // Steps
  const size_t odgi_get_step_count(const odgi::graph_t *graph, const handle_t handle)
  {
    return graph->get_step_count(handle);
  }

  const handle_t odgi_step_get_handle(const graph_t *graph, step_handle_t step)
  {
    return graph->get_handle_of_step(step);
  }

  const path_handle_t odgi_step_get_path(const graph_t *graph, step_handle_t step)
  {
    return graph->get_path(step);
  }

  const step_handle_t odgi_step_path_begin(const graph_t *graph, path_handle_t path)
  {
    return graph->path_begin(path);
  }

  const step_handle_t odgi_step_path_end(const graph_t *graph, path_handle_t path)
  {
    return graph->path_end(path);
  }

  const step_handle_t odgi_step_path_back(const graph_t *graph, path_handle_t path)
  {
    return graph->path_back(path);
  }

  const int64_t odgi_step_path_id(const graph_t *graph, step_handle_t &step_handle) {
    return (reinterpret_cast<const int64_t*>(&step_handle)[0]) >> 1;
  }

  bool odgi_step_is_reverse(const graph_t *graph, step_handle_t &step_handle) {
    return (reinterpret_cast<const int64_t*>(&step_handle)[0]) & 1;
  }

  const int64_t odgi_step_prev_id(const graph_t *graph, step_handle_t &step_handle)
  {
      return (reinterpret_cast<const int64_t*>(&step_handle)[1]);
  }

  const int64_t odgi_step_prev_rank(const graph_t *graph, step_handle_t &step_handle)
  {
      return (reinterpret_cast<const int64_t*>(&step_handle)[2]);
  }

  const int64_t odgi_step_next_id(const graph_t *graph, step_handle_t &step_handle)
  {
      return (reinterpret_cast<const int64_t*>(&step_handle)[3]);
  }

  const int64_t odgi_step_next_rank(const graph_t *graph,step_handle_t &step_handle)
  {
      return (reinterpret_cast<const int64_t*>(&step_handle)[4]);
  }

  const bool odgi_step_eq(const graph_t *graph,
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

  const step_handle_t odgi_path_front_end(const graph_t *graph,
                                     const path_handle_t path)
  {
    return graph->path_front_end(path);
  }

  const step_handle_t odgi_get_next_step(const graph_t *graph,
                                         const step_handle_t step)
  {
    return graph->get_next_step(step);
  }

  const step_handle_t odgi_get_previous_step(const graph_t *graph,
                                         const step_handle_t step)
  {
    return graph->get_previous_step(step);
  }

  const bool odgi_has_edge(const graph_t *graph, const handle_t left, const handle_t right)
  {
    return graph->has_edge(left,right);
  }

  const bool odgi_is_path_front_end(const graph_t *graph, const step_handle_t step)
  {
    return graph->is_path_front_end(step);
  }
  const bool odgi_is_path_end(const graph_t *graph, const step_handle_t step)
  {
    return graph->is_path_end(step);
  }
  const bool odgi_has_next_step(const graph_t *graph, const step_handle_t step)
  {
    return graph->has_next_step(step);
  }
  const bool odgi_has_previous_step(const graph_t *graph, const step_handle_t step)
  {
    return graph->has_previous_step(step);
  }

  const path_handle_t odgi_get_path_handle_of_step(const graph_t *graph,
                                                   const step_handle_t step)
  {
    return graph->get_path_handle_of_step(step);
  }

  void odgi_for_each_step_in_path(const graph_t *graph,
                                  const path_handle_t path,
                                  void (*next) (const step_handle_t step))
  {
    graph->for_each_step_in_path(path,
                                 [&](const step_handle_t step) {
                                   next(step);
                                 });
  };


  const bool odgi_for_each_step_on_handle(const graph_t *graph,
                                          const handle_t handle,
                                          bool (*next) (const step_handle_t step))
  {
    return graph->for_each_step_on_handle(handle,
                                          [&](const step_handle_t step) {
                                            next(step);
                                            return true; });
  }

} // extern "C"
