#pragma once

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <vector>
#include "odgi.hpp"
#include "topological_sort.hpp"
#include "eades_algorithm.hpp"
#include "dagify.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

std::vector<handle_t> cycle_breaking_sort(const graph_t& graph);

}

}
