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

/// Provide depth of our given path ranges to callback, requires the graph to be optimized!
		void for_each_path_range_degree(const PathHandleGraph& graph,
									   const std::vector<path_range_t>& path_ranges,
									   const std::vector<bool>& paths_to_consider,
									   const std::function<void(const path_range_t&, const double&)>& func);

	}

}
