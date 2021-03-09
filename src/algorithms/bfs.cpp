#include "bfs.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;


// breadth first search across handles in the graph
void bfs(
    const HandleGraph& graph,
    const std::function<void(const handle_t&, const uint64_t&, const uint64_t&, const uint64_t&)>& handle_fn,
    const std::function<bool(const handle_t&)>& seen_handle_fn,
    const std::function<bool(const handle_t&, const handle_t&)>& seen_edge_fn,
    const std::function<bool(void)>& break_fn,
    const std::vector<handle_t>& sources,
    const std::vector<handle_t>& sinks,
    bool bidirectional,
    uint64_t step_limit,
    uint64_t bp_limit
    ) {

    ska::flat_hash_set<handle_t> stops; // sinks
    for (auto& handle : sinks) {
        stops.insert(handle);
    }

    std::deque<bfs_state_t> todo; // our traversal front
    uint64_t seed_rank = 0;
    for (auto& handle : sources) {
        todo.push_front({handle, seed_rank++, 0, 0});
    }

    while (!todo.empty()) {
        // get our next handle
        auto& curr = todo.back();
        handle_t handle = curr.handle;
        uint64_t root = curr.root;
        uint64_t curr_length = curr.length;
        uint64_t curr_depth = curr.depth;
        // pop the handle off the back of the queue
        todo.pop_back();
        if (!seen_handle_fn(handle)) {
            // handle the handle
            handle_fn(handle, root, curr_length, curr_depth);
            curr_length += graph.get_length(handle);
            ++curr_depth;
            // check if we've hit our break condition
            if (break_fn()) { return; }
            // check if we should stop here
            if (!stops.count(handle)
                && (!step_limit || curr_depth < step_limit)
                && (!bp_limit || curr_length < bp_limit)
                ) {
                // add the next nodes to our queue
                auto enqueue =
                    [&todo,&handle,&seen_edge_fn,&root,&curr_length,&curr_depth]
                    (const handle_t& next) {
                        if (!seen_edge_fn(handle, next)) {
                            todo.push_back({next, root, curr_length, curr_depth});
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
