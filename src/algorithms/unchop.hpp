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
 * Unchop by gluing abutting handles with just a single edge between them and
 * compatible path steps together.
 */
void unchop(handlegraph::MutablePathDeletableHandleGraph& graph);

//std::vector<std::deque<handle_t>> simple_components(PathHandleGraph* graph, int min_size = 1);

handle_t concat_nodes(handlegraph::MutablePathDeletableHandleGraph& graph, const std::vector<handle_t>& nodes);
//handle_t concat_nodes(handlegraph::MutablePathDeletableHandleGraph* graph, const std::vector<handle_t>& nodes);
    
}
}
