#pragma once

#include <handlegraph/handle_graph.hpp>
#include <vector>
#include <set>
#include <deque>
#include "hash_map.hpp"

namespace odgi {
namespace algorithms {
    
using namespace handlegraph;

void bfs(
    const HandleGraph& graph,
    const std::function<bool(const handle_t&)>& handle_begin_fn,  // called when node orientation is first encountered, return true to continue
    const std::function<bool(void)>& break_fn,                    // called to check if we should stop the DFS; we stop when true is returned.
    //const std::function<void(const edge_t&)>& edge_fn,            // called when an edge is encountered
    const std::vector<handle_t>& sources,                         // start only at these node traversals
    const std::vector<handle_t>& sinks                            // when hitting a sink, don't keep walking
    );

}
}


