#include "coverage.hpp"

namespace odgi {
namespace algorithms {

std::vector<handle_t> find_handles_exceeding_coverage_limits(const MutablePathDeletableHandleGraph& graph, uint64_t min_coverage, uint64_t max_coverage) {
    std::vector<handle_t> handles;
    graph.for_each_handle([&](const handle_t& handle) {
            uint64_t step_count = graph.steps_of_handle(handle).size();
            if (min_coverage && step_count < min_coverage || max_coverage && step_count > max_coverage) {
                handles.push_back(handle);
            }
        });
    return handles;
}

std::vector<edge_t> find_edges_exceeding_coverage_limits(const MutablePathDeletableHandleGraph& graph, uint64_t min_coverage, uint64_t max_coverage) {
    std::vector<edge_t> edges;
    graph.for_each_handle([&](const handle_t& handle) {
            hash_map<handle_t, uint64_t> nexts;
            hash_map<handle_t, uint64_t> prevs;
            graph.for_each_step_on_handle(handle, [&](const step_handle_t& step) {
                    if (graph.has_next_step(step)) {
                        step_handle_t next = graph.get_next_step(step);
                        ++nexts[graph.get_handle_of_step(next)];
                    }
                    if (graph.has_previous_step(step)) {
                        step_handle_t prev = graph.get_previous_step(step);
                        ++prevs[graph.get_handle_of_step(prev)];
                    }
                });
            for (auto& n : nexts) {
                const uint64_t& path_cov = n.second;
                if (min_coverage && path_cov < min_coverage || max_coverage && path_cov > max_coverage) {
                    edges.push_back(std::make_pair(handle, n.first));
                }
            }
            for (auto& p : prevs) {
                const uint64_t& path_cov = p.second;
                if (min_coverage && path_cov < min_coverage || max_coverage && path_cov > max_coverage) {
                    edges.push_back(std::make_pair(p.first, handle));
                }
            }
        });
    std::sort(edges.begin(), edges.end(), [](const edge_t& a, const edge_t& b) {
            return std::make_pair(as_integer(a.first), as_integer(a.second))
                < std::make_pair(as_integer(b.first), as_integer(b.second));
        });
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
}

std::vector<edge_t> keep_mutual_best_edges(const MutablePathDeletableHandleGraph& graph, uint64_t n_best) {
    std::vector<edge_t> edges;
    graph.for_each_handle([&](const handle_t& handle) {
            hash_map<handle_t, uint64_t> nexts;
            hash_map<handle_t, uint64_t> prevs;
            graph.for_each_step_on_handle(handle, [&](const step_handle_t& step) {
                    if (graph.has_next_step(step)) {
                        step_handle_t next = graph.get_next_step(step);
                        ++nexts[graph.get_handle_of_step(next)];
                    }
                    if (graph.has_previous_step(step)) {
                        step_handle_t prev = graph.get_previous_step(step);
                        ++prevs[graph.get_handle_of_step(prev)];
                    }
                });
            std::map<uint64_t, handle_t> nexts_sorted;
            for (auto& n : nexts) {
                nexts_sorted[n.second] = n.first;
            }
            std::map<uint64_t, handle_t> prevs_sorted;
            for (auto& p : prevs) {
                prevs_sorted[p.second] = p.first;
            }
            uint64_t i = 0;
            for (auto& n : nexts_sorted) {
                if (nexts.size() - i++ > n_best) {
                    edges.push_back(std::make_pair(handle, n.second));
                } else {
                    break;
                }
            }
            i = 0;
            for (auto& p : prevs_sorted) {
                if (prevs.size() - i++ > n_best) {
                    edges.push_back(std::make_pair(p.second, handle));
                } else {
                    break;
                }
            }
        });
    std::sort(edges.begin(), edges.end(), [](const edge_t& a, const edge_t& b) {
            return std::make_pair(as_integer(a.first), as_integer(a.second))
                < std::make_pair(as_integer(b.first), as_integer(b.second));
        });
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
}

}
}
