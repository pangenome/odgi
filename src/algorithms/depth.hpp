#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <omp.h>
#include "hash_map.hpp"
#include "position.hpp"
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

/// Find handles with more or less than the given path depth limits
std::vector<handle_t> find_handles_exceeding_depth_limits(const MutablePathDeletableHandleGraph& graph, uint64_t min_depth, uint64_t max_depth);

/// Find edges with more or less than the given path depth limits
std::vector<edge_t> find_edges_exceeding_depth_limits(const MutablePathDeletableHandleGraph& graph, uint64_t min_depth, uint64_t max_depth);

/// Keep the N best edges by path depth inbound and outbound of every node where they are the best for their neighbors
std::vector<edge_t> keep_mutual_best_edges(const MutablePathDeletableHandleGraph& graph, uint64_t n_best);

/// Provide depth of our given path ranges to callback, requires the graph to be optimized!
void for_each_path_range_depth(const PathHandleGraph& graph,
                               const std::vector<path_range_t>& path_ranges,
                               const std::vector<bool>& paths_to_consider,
                               const std::function<void(const path_range_t&, const double&)>& func);

/// Destroy handles with more or less than the given path depth limits
//void bound_depth(MutablePathDeletableHandleGraph& graph, uint64_t min_depth, uint64_t max_depth);

}

}
