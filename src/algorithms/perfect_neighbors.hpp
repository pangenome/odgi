#pragma once

#include <handlegraph/types.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/path_handle_graph.hpp>

namespace odgi {
namespace algorithms {

using namespace handlegraph;

bool nodes_are_perfect_path_neighbors(const PathHandleGraph& graph, handle_t left_handle, handle_t right_handle);

}
}

