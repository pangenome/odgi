#include "path_keep.hpp"

namespace odgi {
namespace algorithms {

void keep_paths(const graph_t& graph, graph_t& into, const ska::flat_hash_set<path_handle_t>& to_keep) {
    graph.for_each_handle([&](const handle_t& h) {
        into.create_handle(graph.get_sequence(h), graph.get_id(h));
    });
    graph.for_each_handle([&](const handle_t& h) {
        graph.follow_edges(h, false, [&](const handle_t& t) {
            into.create_edge(into.get_handle(graph.get_id(h), graph.get_is_reverse(h)),
                             into.get_handle(graph.get_id(t), graph.get_is_reverse(t)));
        });
        auto r = graph.flip(h);
        graph.follow_edges(r, false, [&](const handle_t& t) {
            into.create_edge(into.get_handle(graph.get_id(r), graph.get_is_reverse(r)),
                             into.get_handle(graph.get_id(t), graph.get_is_reverse(t)));
        });
    });
    std::vector<path_handle_t> paths;
    graph.for_each_path_handle([&](const path_handle_t& p) { paths.push_back(p); });
    // ensure we have the same path handle order as in the original graph
    for (auto& path : paths) {
        if (to_keep.count(path)) {
            into.create_path_handle(graph.get_path_name(path));
        }
    }
    // add the steps in parallel
#pragma omp parallel for
    for (auto& path : paths) {
        if (to_keep.count(path)) {
            // add the path as-is
            auto w = into.get_path_handle(graph.get_path_name(path));
            graph.for_each_step_in_path(
                path,
                [&w,&into,&graph](const step_handle_t& s) {
                    auto h = graph.get_handle_of_step(s);
                    auto q = into.get_handle(graph.get_id(h), graph.get_is_reverse(h));
                    into.append_step(w, q);
                });
        }
    }
    into.optimize();
}

}
}
