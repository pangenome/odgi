#pragma once

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
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
    // the second parameter gives the rank of the root among the input set of seeds
    // the third parameter gives the cumulative length from the search root that led to this handle
    // the fourth parameter gives the search depth at this point
    // return true to continue through this node, false to treat it like a sink
    const std::function<void(const handle_t&, const uint64_t&, const uint64_t&, const uint64_t&)>& handle_fn,
    // should we consider this handle and its edges?, returns false if we should continue
    const std::function<bool(const handle_t&)>& seen_handle_fn,
    // the edge we will traverse to get to a new handle, returns false if we should continue across the edge
    const std::function<bool(const handle_t&, const handle_t&)>& seen_edge_fn,
    // called to check if we should stop the BFS; we stop when true is returned.
    const std::function<bool(void)>& break_fn,
    // start only at these node traversals
    const std::vector<handle_t>& sources,
    // when hitting a sink, don't keep walking
    const std::vector<handle_t>& sinks,
    // do we use a bidirectional search
    bool bidirectional = false,
    uint64_t step_limit = 0,
    uint64_t bp_limit = 0
    );

struct bfs_state_t {
    handle_t handle;
    uint64_t root;
    uint64_t length;
    uint64_t depth;
};

}
}


