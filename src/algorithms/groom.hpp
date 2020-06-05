#pragma once

#include <handlegraph/types.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>

#include <vector>

#include "dynamic.hpp"
#include "topological_sort.hpp"
//#include "bfs.cpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

/**
 * Remove spurious inverting links based on a dominant orientation of the graph
 */
void groom(handlegraph::MutablePathDeletableHandleGraph& source,
           handlegraph::MutablePathDeletableHandleGraph& target);

}
}
