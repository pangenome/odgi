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
    // called when a given node orientation is first encountered
    // the second parameter gives the cumulative length from the search root that led to this handle
    // return true to continue through this node, false to treat it like a sink
    const std::function<void(const handle_t&, const uint64_t&)>& handle_fn,
    // have we seen this handle before?
    const std::function<bool(const handle_t&)>& seen_fn,
    // called to check if we should stop the DFS; we stop when true is returned.
    const std::function<bool(void)>& break_fn,
    // start only at these node traversals
    const std::vector<handle_t>& sources,
    // when hitting a sink, don't keep walking
    const std::vector<handle_t>& sinks
    );

}
}


