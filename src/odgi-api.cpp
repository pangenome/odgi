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
//
// To test for memory errors compile with
//
//   cmake -DCMAKE_BUILD_TYPE=Debug -DEXTRA_FLAGS="-fsanitize=address" ..
//   make odgi_ffi -j 16 && ctest . -R ffi --verbose

#include "odgi-api.h"
#include "version.hpp"

using namespace odgi;

const std::string odgi_version() {
  return Version::get_version();
}

size_t odgi_long_long_size() {
  return sizeof(long long)*8;
}

size_t odgi_handle_i_size() {
  return sizeof(handle_i)*8;
}

size_t odgi_step_handle_i_size() {
  return sizeof(step_handle_i)*8;
}

unsigned __int128 odgi_test_uint128() { // FIXME
  // unsigned __int128 ll = (unsigned __int128)0xAB << 64 + 0xCD;
  // return ll;
  return 0;
}


// has python doctest
const ograph_t odgi_load_graph(const char *filen) {
  auto graph = ograph_t(new graph_t());
  std::ifstream in(filen);
  graph->deserialize(in);
  return graph;
}

// has python doctest
void odgi_free_graph(ograph_t graph) {
  // no need to delete a shared ptr!
}

// has python doctest
const size_t odgi_get_node_count(ograph_t graph) {
  return ((as_graph_t(graph)))->get_node_count();
}

// has python doctest
const size_t odgi_max_node_id(const ograph_t graph) {
  return ((as_graph_t(graph)))->max_node_id();
}

const size_t odgi_min_node_id(const ograph_t graph) {
  return ((as_graph_t(graph)))->min_node_id();
}

// has python doctest
const size_t odgi_get_path_count(ograph_t graph) {
  return ((as_graph_t(graph)))->get_path_count();
}

// has python doctest
void odgi_for_each_path_handle(const ograph_t graph,
                               const std::function<void(const path_handle_i)>& next) {
  ((as_graph_t(graph)))->for_each_path_handle([&](const path_handle_t path) {
    next(as_path_handle_i(path));
  });
}

// has python doctest
const bool odgi_for_each_handle(const ograph_t graph,
                                const std::function<bool(const handle_i)>& next)
{
  return ((as_graph_t(graph)))->for_each_handle([&](const handle_t h)
  {

    next(as_handle_i(h));
    return true;
  },
    false); // note that parallel is false
}

/// Loop over all the handles to next/previous (right/left) nodes. Passes
/// them to a callback which returns false to stop iterating and true to
/// continue. Returns true if we finished and false if we stopped early.
// has python doctest
const bool odgi_follow_edges(const ograph_t graph,
                             const handle_i ihandle,
                             bool go_left,
                             const std::function<bool(const handle_i ihandle)>& next)
{
  return ((as_graph_t(graph)))->follow_edges(as_handle(ihandle),go_left,
                             [&](const handle_t handle)
                             {
                               next(as_handle_i(handle));
                               return true;
                             }
                             );
};

const handle_i odgi_edge_first_handle(const ograph_t graph, const edge_t edge_handle)
{
  return as_handle_i((&edge_handle)->first);
}

const handle_i odgi_edge_second_handle(const ograph_t graph, const edge_t edge_handle)
{
  return as_handle_i((&edge_handle)->second);
}

const bool odgi_has_node(const ograph_t graph, nid_t node_id) {
  return ((as_graph_t(graph)))->has_node(node_id);
}

// has python doctest
const std::string odgi_get_sequence(const ograph_t graph, const handle_i ihandle) {
  return ((as_graph_t(graph)))->get_sequence(as_handle(ihandle));
}

// has python doctest
const nid_t odgi_get_id(const ograph_t graph, const handle_i ihandle) {
  return ((as_graph_t(graph)))->get_id(as_handle(ihandle));
}

const auto& odgi_get_node_id = odgi_get_id;

const bool odgi_get_is_reverse(const ograph_t graph, const handle_i ihandle) {
  return ((as_graph_t(graph)))->get_is_reverse(as_handle(ihandle));
}

const size_t odgi_get_length(const ograph_t graph, const handle_i ihandle) {
  return ((as_graph_t(graph)))->get_length(as_handle(ihandle));
}


const bool odgi_has_path(const ograph_t graph, const char *path_name) {
  return ((as_graph_t(graph)))->has_path(path_name);
}

const bool odgi_path_is_empty(const ograph_t graph, const path_handle_i path) {
  return ((as_graph_t(graph)))->is_empty(as_path_handle(path));
}

// has python doctest
const path_handle_i odgi_get_path_handle(const ograph_t graph, const char *path_name) {
  return as_path_handle_i(((as_graph_t(graph)))->get_path_handle(path_name));
}

// Steps
const size_t odgi_get_step_count(const ograph_t graph, const handle_i ihandle)
{
  return ((as_graph_t(graph)))->get_step_count(as_handle(ihandle));
}

// has python doctest
const handle_i odgi_get_handle_of_step(const ograph_t graph, step_handle_i step)
{
  return as_handle_i(((as_graph_t(graph)))->get_handle_of_step(as_step_handle_t(step)));
}

// step -> path
const path_handle_i odgi_get_path(const ograph_t graph, step_handle_i step)
{
  return as_path_handle_i(((as_graph_t(graph)))->get_path(as_step_handle_t(step)));
}

// path -> step
// has python doctest
const step_handle_i odgi_path_begin(const ograph_t graph, path_handle_i path)
{
  return as_step_handle_i(((as_graph_t(graph)))->path_begin(as_path_handle(path)));
}

const step_handle_i odgi_path_end(const ograph_t graph, path_handle_i path)
{
  return as_step_handle_i(((as_graph_t(graph)))->path_end(as_path_handle(path)));
}

const step_handle_i odgi_path_back(const ograph_t graph, path_handle_i ipath)
{
  return as_step_handle_i(((as_graph_t(graph)))->path_back(as_path_handle(ipath)));
}

const int64_t odgi_step_path_id(const ograph_t graph, step_handle_i step_handle) {
  return (reinterpret_cast<const int64_t*>(&step_handle)[0]) >> 1;
}

bool odgi_step_is_reverse(const ograph_t graph, step_handle_i step_handle) {
  return (reinterpret_cast<const int64_t*>(&step_handle)[0]) & 1;
}

const int64_t odgi_step_prev_id(const ograph_t graph, step_handle_i step_handle)
{
    return (reinterpret_cast<const int64_t*>(&step_handle)[1]);
}

const int64_t odgi_step_prev_rank(const ograph_t graph, step_handle_i step_handle)
{
    return (reinterpret_cast<const int64_t*>(&step_handle)[2]);
}

const int64_t odgi_step_next_id(const ograph_t graph, step_handle_i step_handle)
{
    return (reinterpret_cast<const int64_t*>(&step_handle)[3]);
}

const int64_t odgi_step_next_rank(const ograph_t graph,step_handle_i step_handle)
{
    return (reinterpret_cast<const int64_t*>(&step_handle)[4]);
}

const bool odgi_step_eq(const ograph_t graph,
                        step_handle_i step_handle1,
                        step_handle_i step_handle2)
{
  // this code is hardly optimal. FIXME.
  for (uint i=0;i<8;i++) {
    if ((reinterpret_cast<const int64_t*>(&step_handle1)[i]) != (reinterpret_cast<const int64_t*>(&step_handle2)[i])){
      return false;
    }
  }
  return true;
}

const step_handle_i odgi_path_front_end(const ograph_t graph,
                                   const path_handle_i ipath)
{
  return as_step_handle_i(((as_graph_t(graph)))->path_front_end(as_path_handle(ipath)));
}

// has python doctest
const step_handle_i odgi_get_next_step(const ograph_t graph,
                                       const step_handle_i step)
{
  return as_step_handle_i(((as_graph_t(graph)))->get_next_step(as_step_handle_t(step)));
}

const step_handle_i odgi_get_previous_step(const ograph_t graph,
                                       const step_handle_i step)
{
  return as_step_handle_i(((as_graph_t(graph)))->get_previous_step(as_step_handle_t(step)));
}

const bool odgi_has_edge(const ograph_t graph, const handle_i left, const handle_i right)
{
  return ((as_graph_t(graph)))->has_edge(as_handle(left),as_handle(right));
}

const bool odgi_is_path_front_end(const ograph_t graph, const step_handle_i step)
{
  return ((as_graph_t(graph)))->is_path_front_end(as_step_handle_t(step));
}
const bool odgi_is_path_end(const ograph_t graph, const step_handle_i step)
{
  return ((as_graph_t(graph)))->is_path_end(as_step_handle_t(step));
}
const bool odgi_has_next_step(const ograph_t graph, const step_handle_i step)
{
  return ((as_graph_t(graph)))->has_next_step(as_step_handle_t(step));
}
const bool odgi_has_previous_step(const ograph_t graph, const step_handle_i step)
{
  return ((as_graph_t(graph)))->has_previous_step(as_step_handle_t(step));
}

const path_handle_i odgi_get_path_handle_of_step(const ograph_t graph,
                                                 const step_handle_i step)
{
  return as_path_handle_i(((as_graph_t(graph)))->get_path_handle_of_step(as_step_handle_t(step)));
}

void odgi_for_each_step_in_path(const ograph_t graph,
                                const path_handle_i ipath,
                                void (*next) (const step_handle_i step))
{
  ((as_graph_t(graph)))->for_each_step_in_path(as_path_handle(ipath),
                               [&](const step_handle_t step) {
                                 next(as_step_handle_i(step));
                               });
};


const bool odgi_for_each_step_on_handle(const ograph_t graph,
                                        const handle_i ihandle,
                                        bool (*next) (const step_handle_i step))
{
  return ((as_graph_t(graph)))->for_each_step_on_handle(as_handle(ihandle),
                                        [&](const step_handle_t step) {
                                          next(as_step_handle_i(step));
                                          return true; });
}

// has python doctest
const std::string odgi_get_path_name(const ograph_t graph, const path_handle_i ipath) {
  return (as_graph_t(graph))->get_path_name(as_path_handle(ipath));
}
