#pragma once

#include "odgi.hpp"
#include <iostream>
#include <omp.h>
#include "position.hpp"
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/mutable_path_handle_graph.hpp>
#include "hash_map.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

void keep_paths(const graph_t& graph, graph_t& into, const ska::flat_hash_set<path_handle_t>& to_keep);

}
}
