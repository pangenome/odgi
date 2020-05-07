#pragma once

#include <handlegraph/types.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>

#include <vector>

#include "unchop.hpp"
#include "simplify_siblings.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

/**
 * Normalize a graph, performing up to the given number of iterations.
 * Simplifies siblings and unchops runs of nodes, in a loop.
 */
void normalize(handlegraph::MutablePathDeletableHandleGraph& graph, int max_iter = 1, bool debug = false);
    
}
}
