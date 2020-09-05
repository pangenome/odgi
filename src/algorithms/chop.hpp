#pragma once

#include <handlegraph/types.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>

#include <vector>

#include "simple_components.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

/**
 * Cut nodes to be less than the given max node length.
 */
void chop(handlegraph::MutablePathDeletableHandleGraph& graph, const uint64_t& max_node_length,
          const uint64_t& nthreads, const bool& show_info);
    
}
}
