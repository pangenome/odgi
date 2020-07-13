#pragma once

#include <handlegraph/types.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>

#include <vector>

#include "hash_map.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

/**
 * Merge the given ranges of bases on the given handles together, rewriting paths.
 * Sequences must match. Handles to a single node may occur no more than once.
 */
void merge(handlegraph::MutablePathDeletableHandleGraph& graph, const std::vector<std::pair<handle_t, size_t>>& start, size_t length);
    
}
}

