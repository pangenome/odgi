#include "path_keep.hpp"

namespace odgi {
namespace algorithms {

void keep_paths(const graph_t& graph, graph_t& into, const ska::flat_hash_set<path_handle_t>& to_keep) {
    // for each path, find its average orientation
    graph.for_each_handle([&](const handle_t& h) {
        into.create_handle(graph.get_sequence(h), graph.get_id(h));
    });
    std::vector<path_handle_t> paths;
    graph.for_each_path_handle([&](const path_handle_t& p) { paths.push_back(p); });
#pragma omp parallel for
    for (auto& path : paths) {
        if (to_keep.count(path)) {
            // add the path as-is
            auto w = into.create_path_handle(graph.get_path_name(path));
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
