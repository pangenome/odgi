#include "bfs.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;


// breadth first search across handles in the graph
void bfs(
    const HandleGraph& graph,
    // called when a given node orientation is first encountered
    // the second parameter gives the cumulative length from the search root that led to this handle
    const std::function<void(const handle_t&, const uint64_t&)>& handle_fn,
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

    struct handle_length_t {
        handle_t handle;
        uint64_t length;
    };
    std::deque<handle_length_t> todo; // our traversal front
    for (auto& handle : sources) {
        todo.push_front({handle, 0});
    }

    while (!todo.empty()) {
        // get our next handle
        handle_t handle = todo.back().handle;
        uint64_t curr_length = todo.back().length + graph.get_length(handle);
        // pop the handle off the back of the queue
        todo.pop_back();
        if (!seen_fn(handle)) {
            // handle the handle
            handle_fn(handle, curr_length);
            // check if we've hit our break condition
            if (break_fn()) { return; }
            // check if we should stop here
            if (!stops.count(handle)) {
                // add the next nodes to our queue
                auto enqueue
                    = [&todo,&seen_fn,&curr_length](const handle_t& next) {
                          if (!seen_fn(next)) {
                              todo.push_back({next, curr_length});
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
