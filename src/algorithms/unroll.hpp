#pragma once

#include <handlegraph/types.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
#include <vector>
#include <unordered_set>
#include <list>
#include <set>
#include <iostream>
#include <sstream>
#include <atomic>
#include "hash_map.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

/**
 * Apply linear unrolling to linearize looping components in the graph
 */
void unroll(
    const handlegraph::PathHandleGraph& input,
    handlegraph::MutablePathDeletableHandleGraph& output);

}
}
