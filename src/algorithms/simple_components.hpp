#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include "threads.hpp"
#include <handlegraph/types.hpp>
#include <handlegraph/iteratee.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/mutable_handle_graph.hpp>
#include <handlegraph/mutable_path_handle_graph.hpp>
#include <handlegraph/mutable_path_mutable_handle_graph.hpp>
#include <handlegraph/deletable_handle_graph.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>

namespace odgi {

namespace algorithms {

using namespace handlegraph;

std::vector<std::vector<handle_t>> simple_components(const HandleGraph& graph, uint64_t min_size);

}

}
