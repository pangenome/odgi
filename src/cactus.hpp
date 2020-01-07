#pragma once

/** \file
 * cactus.hpp: Wrapper utility functions for for pinchesAndCacti
 */

#include <vector>
#include <map>
#include <utility>
#include <string>
#include <functional>
#include <queue>

#include "odgi.hpp"
#include "hash_map.hpp"

#include "topological_sort.hpp"
#include "weakly_connected_components.hpp"
#include "strongly_connected_components.hpp"
#include "find_shortest_paths.hpp"
#include "dfs.hpp"

//#include "types.hpp"
//#include "utility.hpp"
//#include "nodeside.hpp"
#include <handlegraph/path_handle_graph.hpp>

extern "C" {
#include "sonLib.h"
#include "stCactusGraphs.h"
}

namespace odgi {

using namespace handlegraph;

// We use this when talking to Cactus.
struct CactusSide {
    int64_t node;
    bool is_end;
};

// Convert VG to Cactus Graph. Takes a list of path names to use to find
// telomeres if present in a connected component.
// Notes:
//  - returned cactus graph needs to be freed by stCactusGraph_destruct
//  - returns a Cactus graph, and a list of stCactusEdgeEnd* telomeres, in std::pairs of adjacent items.
std::pair<stCactusGraph*, stList*> handle_graph_to_cactus(PathHandleGraph& graph, const ska::flat_hash_set<std::string>& hint_paths);

// Convert back from Cactus to VG
// (to, for example, display using vg view)
// todo: also provide mapping info to get nodes embedded in cactus components
//VG cactus_to_vg(stCactusGraph* cactus_graph);

// Convert vg into vg formatted cactus representation
// Input graph must be sorted!
//VG cactusify(VG& graph);

}

