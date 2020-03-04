#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "threads.hpp"
#include "hash_map.hpp"
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

/// Find edges with more or less than the given path coverage limits
std::vector<edge_t> find_edges_exceeding_coverage_limits(const MutablePathDeletableHandleGraph& graph, uint64_t min_coverage, uint64_t max_coverage);

/// Keep the N best edges by path coverage inbound and outbound of every node where they are the best for their neighbors
std::vector<edge_t> keep_mutual_best_edges(const MutablePathDeletableHandleGraph& graph, uint64_t n_best);

/// Destroy handles with more or less than the given path coverage limits
//void bound_coverage(MutablePathDeletableHandleGraph& graph, uint64_t min_coverage, uint64_t max_coverage);

}

}
