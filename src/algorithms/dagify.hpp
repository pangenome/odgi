#pragma once

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/mutable_handle_graph.hpp>
#include "hash_map.hpp"
#include "subgraph.hpp"
#include "shortest_cycle.hpp"
#include "eades_algorithm.hpp"
#include "strongly_connected_components.hpp"
#include "is_single_stranded.hpp"
#include "topological_sort.hpp"
//#include <unordered_map>

namespace odgi {

namespace algorithms {

// Fill an empty MutableHandleGraph with a copy of graph where nodes and edges have
// been duplicated in such a way as to eliminate cycles while preserving all paths
// up to a given minimum length. Input HandleGraph must have a single stranded orientation.
// Consider checking this property with has_single_stranded_orientation() before using.
// Returns a mapping from the node IDs of into to the node IDs in graph.
ska::flat_hash_map<handlegraph::nid_t, handlegraph::nid_t> dagify(const HandleGraph* graph, MutableHandleGraph* into,
                                                                  size_t min_preserved_path_length);
}
}
