#include "bfs.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;


// breadth first search across handles in the graph
void bfs(
    const HandleGraph& graph,
    const std::function<bool(const handle_t&)>& handle_begin_fn,  // called when node orientation is first encountered, return true to continue
    const std::function<bool(void)>& break_fn,                    // called to check if we should stop the DFS; we stop when true is returned.
    //const std::function<void(const edge_t&)>& edge_fn,            // called when an edge is encountered
    const std::vector<handle_t>& sources,                         // start only at these node traversals
    const std::vector<handle_t>& sinks                            // when hitting a sink, don't keep walking
    ) {

    ska::flat_hash_set<handle_t> seen;  // handles we've seen

    ska::flat_hash_set<handle_t> stops; // sinks
    for (auto& handle : sinks) {
        stops.insert(handle);
    }

    std::deque<handle_t> todo; // our traversal front
    for (auto& handle : sources) {
        todo.push_front(handle);
    }

    while (!todo.empty()) {
        // get our next handle
        handle_t handle = todo.back();
        // pop the handle off the back of the queue
        todo.pop_back();
        if (!seen.count(handle)) {
            // record that we've visited the node
            seen.insert(handle);
            // handle the handle
            if (handle_begin_fn(handle)) {
                // check if we've hit our break condition
                if (break_fn()) { return; }
                // check if we should stop here
                if (!stops.count(handle)) {
                    // add the next nodes to our queue
                    graph.follow_edges(
                        handle, false,
                        [&todo,&seen](const handle_t& next) {
                            if (!seen.count(next)) {
                                todo.push_back(next);
                            }
                        });
                }
            }
        }
    }
}

}
}
