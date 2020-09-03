#pragma once

/**
 * \file weakly_connected_components.hpp
 *
 * Defines an algorithm for finding weakly connected components in a graph.
 */

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include "hash_map.hpp"
#include <vector>
#include <algorithm>

namespace odgi {
namespace algorithms {

using namespace handlegraph;

/// Returns sets of IDs defining components that are connected by any series
/// of nodes and edges, even if it is not a valid bidirected walk. TODO: It
/// might make sense to have a handle-returning version, but the consumers of
/// weakly connected components right now want IDs, and membership in a weakly
/// connected component is orientation-independent.
std::vector<ska::flat_hash_set<handlegraph::nid_t>> weakly_connected_components(const HandleGraph* graph);

/// Returns a vector of handles, one for each component, which can be easier to use in some cases
std::vector<std::vector<handlegraph::handle_t>> weakly_connected_component_vectors(const HandleGraph* graph);

/// Return pairs of weakly connected component ID sets and the handles that are
/// their tips, oriented inward. If a node is both a head and a tail, it will
/// appear in tips in both orientations.
std::vector<std::pair<ska::flat_hash_set<handlegraph::nid_t>, std::vector<handle_t>>> weakly_connected_components_with_tips(const HandleGraph* graph);

}
}
