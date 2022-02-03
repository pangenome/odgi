#include "heaps.hpp"

namespace odgi {

namespace algorithms {

void for_each_heap_permutation(const PathHandleGraph& graph,
                               const std::vector<std::vector<path_handle_t>>& path_groups,
                               uint64_t n_permutations,
                               const std::function<void(const std::vector<uint64_t>&, uint64_t)>& func) {
    //std::vector<std::vector<path_handle_t>>
    auto get_permutation = [&](void) {
        std::random_device rd;
        std::default_random_engine rng(rd());
        std::vector<uint64_t> order; order.reserve(path_groups.size());
        for (uint64_t i = 0; i < path_groups.size(); ++i) {
            order.push_back(i);
        }
        std::shuffle(order.begin(), order.end(), rng);
        return order;
    };

    uint64_t graph_total_bp = 0;
    uint64_t idx = 0;
    graph.for_each_handle([&](const handle_t& h) {
        graph_total_bp += graph.get_length(h);
        if (graph.get_id(h) != ++idx) {
            std::cerr << "[odgi::algorithms::heaps] error: the node IDs are not compacted. Please run 'odgi sort' using -O, --optimize to optimize the graph." << std::endl;
            std::abort();
        }
    });

#pragma omp parallel for
    for (uint64_t i = 0; i < n_permutations; ++i) {
        auto permutation = get_permutation();
        // XXX the graph must be node-id compacted for this trick!
        std::vector<bool> seen_nodes(graph.get_node_count(), 0);
        uint64_t seen_bp = 0;
        std::vector<uint64_t> vals;
        for (auto& j : permutation) {
            auto& paths = path_groups[j];
            for (auto& path : paths) {
                graph.for_each_step_in_path(
                    path,
                    [&](const step_handle_t& step) {
                        auto h = graph.get_handle_of_step(step);
                        uint64_t rank = graph.get_id(h)-1; // assumes compaction!
                        if (!seen_nodes[rank]) {
                            seen_nodes[rank] = 1;
                            seen_bp += graph.get_length(h);
                        }
                    });
            }
            vals.push_back(seen_bp);
        }
        func(vals, i);
    }
}

}

}
