#pragma once

#include <unordered_set>
#include <handlegraph/handle_graph.hpp>
#include "hash_map.hpp"
#include "dfs.hpp"

namespace odgi {
namespace algorithms {

using namespace std;
using namespace handlegraph;

/// Find all of the nodes with no edges on their left sides.
vector<ska::flat_hash_set<handlegraph::nid_t>> strongly_connected_components(const HandleGraph* g);
    
}
}
