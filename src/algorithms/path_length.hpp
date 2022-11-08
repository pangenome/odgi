#pragma once

//#include "odgi.hpp"
#include <iostream>
#include <omp.h>
#include "position.hpp"
#include <handlegraph/path_handle_graph.hpp>
#include "hash_map.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

ska::flat_hash_map<path_handle_t, uint64_t> get_path_length(const PathHandleGraph& graph);

}
}
