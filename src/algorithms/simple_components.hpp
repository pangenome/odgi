#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include "threads.hpp"
#include <handlegraph/types.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/path_handle_graph.hpp>

#include "perfect_neighbors.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

std::vector<std::vector<handle_t>> simple_components(const PathHandleGraph& graph, uint64_t min_size);

}

}
