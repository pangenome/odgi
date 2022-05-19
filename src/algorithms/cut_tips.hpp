#pragma once

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/deletable_handle_graph.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <vector>
#include "odgi.hpp"
#include "bfs.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

std::vector<handle_t> tip_handles(
    const HandleGraph& graph);

uint64_t cut_tips(
    DeletableHandleGraph& graph);

uint64_t cut_tips(
    MutablePathDeletableHandleGraph& graph,
    uint64_t min_depth = 0);

}

}
