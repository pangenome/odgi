#include "cut_tips.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

// find and remove isolated nodes which have only a single path on them
std::vector<std::pair<path_handle_t, handle_t>> isolated_path_handles(
    const MutablePathDeletableHandleGraph& graph) {
    std::vector<std::pair<path_handle_t, handle_t>> isolated;
    graph.for_each_handle(
        [&](const handle_t& handle) {
            uint64_t step_count = 0;
            step_handle_t single_step;
            size_t degree = graph.get_degree(handle, false) + graph.get_degree(handle, true);
            if (degree == 0) {
                graph.for_each_step_on_handle(
                    handle,
                    [&](const step_handle_t& step) {
                        single_step = step;
                        ++step_count;
                    });
                if (step_count == 1) {
                    // is the single step the only step in that path?
                    path_handle_t path = graph.get_path_handle_of_step(single_step);
                    if (single_step == graph.path_begin(path)
                        && single_step == graph.path_back(path)) {
                        isolated.push_back(std::make_pair(path, handle));
                    }
                }
            }
        });
    return isolated;
}

uint64_t remove_isolated_paths(
    MutablePathDeletableHandleGraph& graph) {
    auto isolated = isolated_path_handles(graph);
    for (auto& p : isolated) {
        auto& path = p.first;
        auto& handle = p.second;
        graph.destroy_path(path);
        graph.destroy_handle(handle);
    }
    return isolated.size();
}

}

}
