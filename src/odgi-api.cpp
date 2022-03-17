// odgi-api exports the C-API for FFIs
//
// Copyright (c) 2022 Erik Garrison and Pjotr Prins
//
// ODGI is published under the MIT license

#include "odgi-api.h"
#include "version.hpp"


extern "C" {

  using namespace odgi;

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

  const char *odgi_version() {
    return Version::get_version().c_str();
  }

  void odgi_for_each_path_handle(const graph_t *graph, void (*next) (const path_handle_t path)) {
    graph->for_each_path_handle([&](const path_handle_t& path) {
      next(path);
    });
  }

  const bool odgi_for_each_handle(const graph_t *graph,
                                  // const std::function<bool(const handlegraph::handle_t&)>& iteratee)
                                  bool (*next) (const handlegraph::handle_t handle))
  {
    return graph->for_each_handle([&](const handlegraph::handle_t h)
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
                               const handlegraph::handle_t handle,
                               bool go_left,
                               bool (*next) (const handlegraph::handle_t handle))
  {
    return graph->follow_edges(handle,go_left,
                               [&](const handlegraph::handle_t handle)
                               {
                                 next(handle);
                                 return true;
                               }
                               );
  };

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

  const char *odgi_get_path_name(const graph_t *graph, const path_handle_t path) {
    return graph->get_path_name(path).c_str();
  }

  const bool odgi_has_path(const graph_t *graph, const char *path_name) {
    return graph->has_path(path_name);
  }

  const path_handle_t odgi_get_path_handle(const graph_t *graph, const char *path_name) {
    return graph->get_path_handle(path_name);
  }
}
