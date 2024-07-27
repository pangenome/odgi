#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <random>
#include <omp.h>
#include "hash_map.hpp"
#include "position.hpp"
#include "odgi.hpp"
#include <handlegraph/types.hpp>
#include <handlegraph/iteratee.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>

namespace odgi {

namespace algorithms {

using namespace handlegraph;

/// Flip the orientation of paths in the graph so that they tend to traverse
/// the graph in the forward orientation.
void flip_paths(graph_t& graph,
                graph_t& into,
                const std::vector<path_handle_t>& no_flips,
                const std::vector<path_handle_t>& ref_flips);

}

}
