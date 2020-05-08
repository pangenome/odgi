#include "simple_components.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

// the set of components that could be merged into single nodes without
// changing the path space of the graph
std::vector<std::vector<handle_t>> simple_components(const PathHandleGraph& graph, uint64_t min_size) {
    std::unordered_set<handle_t> seen;
    std::set<std::vector<uint64_t>> components;
    graph.for_each_handle([&](const handle_t& handle) {
            if (!seen.count(handle)) {
                // go left and right through each as far as we have only single edges connecting us
                // to nodes that have only single edges coming in or out
                // that go to other nodes
                std::vector<handle_t> right_linear_component;
                std::unordered_set<handle_t> todo;
                // check the inbound edge count of each node in the direction we are walking
                // if it's > 1 then don't pass through it
                bool stop = false;
                todo.insert(handle);
                while (!stop && todo.size() == 1) {
                    handle_t curr = *todo.begin();
                    if (seen.count(curr)) { stop = true; break; }
                    seen.insert(curr);
                    todo.clear();
                    if (curr != handle) right_linear_component.push_back(curr);
                    graph.follow_edges(curr, false, [&](const handle_t& next) {
                            uint64_t left_edge_count = 0;
                            graph.follow_edges(next, true, [&](const handle_t& h) { ++left_edge_count; });
                            if (graph.get_is_reverse(handle) == graph.get_is_reverse(next)
                                && left_edge_count == 1
                                && nodes_are_perfect_path_neighbors(graph, curr, next)) {
                                todo.insert(next);
                            } else {
                                stop = true;
                            }
                        });
                }
                stop = false;
                todo.clear();
                todo.insert(handle);
                std::vector<handle_t> left_linear_component;
                while (!stop && todo.size() == 1) {
                    handle_t curr = *todo.begin();
                    if (seen.count(curr)) { stop = true; break; }
                    seen.insert(curr);
                    todo.clear();
                    if (curr != handle) left_linear_component.push_back(curr);
                    graph.follow_edges(curr, true, [&](const handle_t& prev) {
                            uint64_t right_edge_count = 0;
                            graph.follow_edges(prev, false, [&](const handle_t& h) { ++right_edge_count; });
                            if (graph.get_is_reverse(handle) == graph.get_is_reverse(prev)
                                && right_edge_count == 1
                                && nodes_are_perfect_path_neighbors(graph, prev, curr)) {
                                todo.insert(prev);
                            } else {
                                stop = true;
                            }
                        });
                }
                std::vector<uint64_t> linear_component;
                linear_component.reserve(left_linear_component.size() + 1 + right_linear_component.size());
                for (auto i = left_linear_component.rbegin(); i != left_linear_component.rend(); ++i) {
                    linear_component.push_back(as_integer(*i));
                }
                linear_component.push_back(as_integer(handle));
                for (auto i = right_linear_component.begin(); i != right_linear_component.end(); ++i) {
                    linear_component.push_back(as_integer(*i));
                }
                if (linear_component.size() >= min_size) {
                    components.insert(linear_component);
                }
            }
        });
    std::vector<std::vector<handle_t>> handle_components;
    for (auto& v : components) {
        handle_components.emplace_back();
        auto& n = handle_components.back();
        for (auto& i : v) {
            n.push_back(as_handle(i));
        }
    }
    return handle_components;
}

}

}
