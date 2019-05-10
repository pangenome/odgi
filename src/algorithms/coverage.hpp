#pragma once

#include <iostream>
#include <vector>
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

using namespace handlegraph;

namespace algorithms {

/// Find handles with more or less than the given path coverage limits
std::vector<handle_t> find_handles_exceeding_coverage_limits(const MutablePathDeletableHandleGraph& graph, uint64_t min_coverage, uint64_t max_coverage);

/// Destroy handles with more or less than the given path coverage limits
//void bound_coverage(MutablePathDeletableHandleGraph& graph, uint64_t min_coverage, uint64_t max_coverage);

}

}
