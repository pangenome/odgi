#include "heaps.hpp"

namespace odgi {

namespace algorithms {

void for_each_heap_permutation(const PathHandleGraph& graph,
                               const std::vector<std::vector<path_handle_t>>& path_groups,
                               const ska::flat_hash_map<path_handle_t, std::vector<interval_t>>& path_intervals,
                               uint64_t n_permutations,
                               uint64_t min_node_depth,
                               const std::function<void(const std::vector<uint64_t>&, uint64_t)>& func) {
    //const std::function<bool(const path_handle_t&, _t)>& in_range) {
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

    // which nodes are traversed by our target paths?
    atomicbitvector::atomic_bv_t target_nodes(graph.get_node_count() + 1);

    if (path_intervals.size() == 0) {
        // keep everything if we aren't given intervals to guide subset
        for (uint64_t i = 0; i < graph.get_node_count(); ++i) {
            if (min_node_depth == 0
                || (graph.get_step_count(graph.get_handle(i+1))
                    >= min_node_depth)) {
                target_nodes.set(i, true);
            }
        }
    } else {
        std::vector<path_handle_t> paths;
        graph.for_each_path_handle([&](const path_handle_t& path) {
            paths.push_back(path);
        });
#pragma omp parallel for
        for (auto& path : paths) {
            if (path_intervals.find(path) != path_intervals.end()) {
                auto& intervals = path_intervals.find(path)->second;
                auto ival = intervals.begin();
                std::set<uint64_t> interval_ends;
                uint64_t pos = 0;
                graph.for_each_step_in_path(
                    path,
                    [&](const step_handle_t& step) {
                        // remove intervals that ended before this node
                        while (interval_ends.size()
                               && *interval_ends.begin() < pos) {
                            interval_ends.erase(interval_ends.begin());
                        }
                        auto h = graph.get_handle_of_step(step);
                        auto len = graph.get_length(h);
                        // add intervals that start on this node
                        while (ival != intervals.end()
                               && ival->first >= pos
                               && ival->first < pos + len) {
                            interval_ends.insert(ival->second);
                            ++ival;
                        }
                        uint64_t rank = graph.get_id(h)-1; // assumes compaction!
                        if (interval_ends.size()
                            && (min_node_depth == 0
                                || graph.get_step_count(h) >= min_node_depth)) {
                            target_nodes.set(rank, true);
                        }
                        pos += len;
                    });
            }
        }
    }

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
                uint64_t pos = 0;
                graph.for_each_step_in_path(
                    path,
                    [&](const step_handle_t& step) {
                        auto h = graph.get_handle_of_step(step);
                        uint64_t rank = graph.get_id(h)-1; // assumes compaction!
                        if (target_nodes[rank] && !seen_nodes[rank]) {
                            seen_nodes[rank] = 1;
                            seen_bp += graph.get_length(h);
                        }
                        pos += graph.get_length(h);
                    });
            }
            vals.push_back(seen_bp);
        }
        func(vals, i);
    }

}

}

}
