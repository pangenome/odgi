#pragma once

//#include <handlegraph/handle_graph.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
//#include <handlegraph/path_handle_graph.hpp>
//#include <handlegraph/deletable_handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <vector>
#include "odgi.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

std::vector<std::pair<path_handle_t, handle_t>> isolated_path_handles(
    const MutablePathDeletableHandleGraph& graph);

uint64_t remove_isolated_paths(
    MutablePathDeletableHandleGraph& graph);

}

}
