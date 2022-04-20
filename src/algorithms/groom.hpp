#pragma once

#include <handlegraph/types.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>

#include <vector>

#include "dynamic.hpp"
#include "topological_sort.hpp"
#include "progress.hpp"
#include "dfs.hpp"

namespace odgi {
    namespace algorithms {

        using namespace handlegraph;

/**
 * Remove spurious inverting links based on a dominant orientation of the graph
 */
        std::vector<handle_t>
        groom(const handlegraph::MutablePathDeletableHandleGraph &graph,
              bool progress_reporting,
			  const std::vector<handlegraph::path_handle_t> target_paths,
              bool use_bfs = true);

    }
}
