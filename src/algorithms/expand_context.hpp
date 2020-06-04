#pragma once

/**
 * \file expand_context.hpp
 *
 * Defines algorithm for adding graph material from the context around a subgraph
 * into the subgraph
 */

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/mutable_handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/mutable_path_mutable_handle_graph.hpp>
#include <handlegraph/util.hpp>

#include "structures/rank_pairing_heap.hpp"

#include <queue>
#include <unordered_set>
#include "hash_map.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

// basic method to query regions of the graph. subgraph must have the same
// ID and sequence space as source.
// use_steps flag toggles whether dist refers to steps or length in base pairs
// note: neighboring nodes are considered to be length 0 away from each other,
// so to do a null expansion by length, use dist < 0.
void expand_context(const HandleGraph* source, MutableHandleGraph* subgraph,
                    int64_t dist, bool use_steps = true, bool expand_forward = true,
                    bool expand_backward = true);
    
// same as above but with addition of paths
// if the subgraph contains disconnected regions of a path, it will attempt to
// add both regions as separate paths with a disambiguated name
void expand_context_with_paths(const PathHandleGraph* source,
                               MutablePathMutableHandleGraph* subgraph,
                               int64_t dist, bool use_steps = true, bool expand_forward = true,
                               bool expand_backward = true);

}
}
