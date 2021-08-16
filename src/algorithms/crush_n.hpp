#pragma once

#include <handlegraph/types.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
#include <vector>
#include "odgi.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

/**
 * Replace runs of Ns at the start and end of nodes with a single N.
 */
void crush_n(odgi::graph_t& graph);
    
}
}
