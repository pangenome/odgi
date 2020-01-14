#include "bfs.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;


// breadth first search across handles in the graph
void bfs(
    const HandleGraph& graph,
    // called when a given node orientation is first encountered
    // the second parameter gives the rank of the root among the input set of seeds
    // the third parameter gives the cumulative length from the search root that led to this handle
    // return true to continue through this node, false to treat it like a sink
    const std::function<void(const handle_t&, const uint64_t&, const uint64_t&)>& handle_fn,
    // have we seen this handle before?
    const std::function<bool(const handle_t&)>& seen_fn,
    // called to check if we should stop the DFS; we stop when true is returned.
    const std::function<bool(void)>& break_fn,
    // start only at these node traversals
    const std::vector<handle_t>& sources,
    // when hitting a sink, don't keep walking
    const std::vector<handle_t>& sinks,
    // do we use a bidirectional search
    bool bidirectional
    ) {

    ska::flat_hash_set<handle_t> stops; // sinks
    for (auto& handle : sinks) {
        stops.insert(handle);
    }

    std::deque<bfs_state_t> todo; // our traversal front
    uint64_t seed_rank = 0;
    for (auto& handle : sources) {
        todo.push_front({handle, seed_rank++, 0});
    }

    while (!todo.empty()) {
        // get our next handle
        auto& curr = todo.back();
        handle_t handle = curr.handle;
        uint64_t root = curr.root;
        uint64_t curr_length = curr.length + graph.get_length(handle);
        // pop the handle off the back of the queue
        todo.pop_back();
        if (!seen_fn(handle)) {
            // handle the handle
            handle_fn(handle, root, curr_length);
            // check if we've hit our break condition
            if (break_fn()) { return; }
            // check if we should stop here
            if (!stops.count(handle)) {
                // add the next nodes to our queue
                auto enqueue =
                    [&todo,&seen_fn,&root,&curr_length]
                    (const handle_t& next) {
                          if (!seen_fn(next)) {
                              todo.push_back({next, root, curr_length});
                          }
                      };
                graph.follow_edges(handle, false, enqueue);
                if (bidirectional) {
                    graph.follow_edges(handle, true, enqueue);
                }
            }
        }
    }
}

}
}
