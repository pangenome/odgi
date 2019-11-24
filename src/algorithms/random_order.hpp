#pragma once

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/mutable_handle_graph.hpp>
#include <algorithm>
#include <random>

namespace odgi {

namespace algorithms {

using namespace handlegraph;

// provide a randomized order for the graph
std::vector<handle_t> random_order(const HandleGraph& graph);

}
}
