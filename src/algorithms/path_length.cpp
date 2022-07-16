#include "path_length.hpp"

namespace odgi {
namespace algorithms {

ska::flat_hash_map<path_handle_t, uint64_t> get_path_length(const PathHandleGraph& graph) {
    ska::flat_hash_map<path_handle_t, uint64_t> path_length;
    std::vector<path_handle_t> paths;
    graph.for_each_path_handle([&](const path_handle_t& path) { paths.push_back(path); });
#pragma omp parallel for
    for (auto& path : paths) {
        uint64_t length = 0;
        graph.for_each_step_in_path(
            path,
            [&](const step_handle_t& step) {
                length += graph.get_length(
                    graph.get_handle_of_step(step));
            });
#pragma omp critical (path_length)
        path_length[path] = length;
    }
    return path_length;
}

}
}
