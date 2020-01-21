#pragma once

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <vector>
#include "odgi.hpp"
#include "bfs.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

std::vector<edge_t> edges_inducing_cycles(
    const HandleGraph& graph,
    const uint64_t& max_cycle_size,
    const uint64_t& max_search_bp
    );

}

}
