#pragma once

#include <handlegraph/types.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
#include <map>
#include "merge.hpp"
#include "progress.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

/**
 * Simplify siblings in the given graph.
 *
 * When one base has two successors with the same base value, and those
 * successors have the same set of predecessors, the successors will be merged.
 *
 * Performs only a subset of the possible merges. Can only merge in from one
 * side of a given node in a single invocation. Returns true if it made
 * progress and there may be more merging to do.
 *
 * Preserves paths.
 */
bool simplify_siblings(handlegraph::MutablePathDeletableHandleGraph& graph,
                       const std::string &progress_message = "");
    
}
}
