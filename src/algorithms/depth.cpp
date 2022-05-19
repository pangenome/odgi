#include "depth.hpp"

namespace odgi {
namespace algorithms {

std::vector<handle_t> find_handles_exceeding_depth_limits(const MutablePathDeletableHandleGraph& graph, uint64_t min_depth, uint64_t max_depth) {
    std::vector<handle_t> handles;
    graph.for_each_handle([&](const handle_t& handle) {
            uint64_t step_count = graph.steps_of_handle(handle).size();
            if (min_depth && step_count < min_depth || max_depth && step_count > max_depth) {
                handles.push_back(handle);
            }
        });
    return handles;
}

std::vector<edge_t> find_edges_exceeding_depth_limits(const MutablePathDeletableHandleGraph& graph, uint64_t min_depth, uint64_t max_depth) {
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
                if (min_depth && path_cov < min_depth || max_depth && path_cov > max_depth) {
                    edges.push_back(std::make_pair(handle, n.first));
                }
            }
            for (auto& p : prevs) {
                const uint64_t& path_cov = p.second;
                if (min_depth && path_cov < min_depth || max_depth && path_cov > max_depth) {
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

void for_each_path_range_depth(const PathHandleGraph& graph,
                               const std::vector<path_range_t>& _path_ranges,
                               const std::vector<bool>& paths_to_consider,
                               const std::function<void(const path_range_t&, const double&)>& func) {
	const uint64_t shift = graph.min_node_id();
	if (graph.max_node_id() - shift >= graph.get_node_count()){
		std::cerr << "[depth::for_each_path_range_depth] error: the node IDs are not compacted. Please run 'odgi sort' using -O, --optimize to optimize the graph." << std::endl;
		exit(1);
	}
    // are we subsetting?
    const bool subset_paths = !paths_to_consider.empty();
    // sort path ranges
    std::vector<const path_range_t*> path_ranges;
    for (auto& p : _path_ranges) {
        path_ranges.push_back(&p);
    }
    std::sort(path_ranges.begin(),
              path_ranges.end(),
              [&](const path_range_t* a,
                  const path_range_t* b) {
                  return std::tie(as_integer(a->begin.path), a->begin.offset, a->end.offset, a->begin.is_rev)
                      < std::tie(as_integer(b->begin.path), b->begin.offset, b->end.offset, b->begin.is_rev);
              });
    // order our path ranges
    std::vector<std::pair<std::vector<const path_range_t*>::const_iterator,
                          std::vector<const path_range_t*>::const_iterator>> paths_todo;
    auto p = path_ranges.begin();
    while (p != path_ranges.end()) {
        auto n = p; ++n;
        while (n != path_ranges.end()
               && (*n)->begin.path == (*p)->begin.path) {
            ++n;
        }
        paths_todo.push_back(std::make_pair(p, n));
        p = n;
    }
    // precompute depths for all handles in parallel
    std::vector<uint64_t> depths(graph.get_node_count() + 1);
    graph.for_each_handle(
        [&](const handle_t& h) {
            auto id = graph.get_id(h);
            auto& d = depths[id - shift];
            if (subset_paths) {
                graph.for_each_step_on_handle(
                    h,
                    [&](const step_handle_t &s) {
                        if (paths_to_consider[as_integer(graph.get_path_handle_of_step(s))]) {
                            ++d;
                        }
                    });
            } else {
                d = graph.get_step_count(h);
            }
        }, true);
    // dip into the path ranges, for each
#pragma omp parallel for schedule(dynamic,1)
    for (auto& p : paths_todo) {
        auto& path = (*p.first)->begin.path;
        std::vector<std::pair<const path_range_t*, uint64_t>> active_ranges;
        auto range_itr = p.first;
        auto ranges_end = p.second;
        uint64_t offset = 0;
        graph.for_each_step_in_path(
            path,
            [&](const step_handle_t& step) {
                // understand our context
                handle_t handle = graph.get_handle_of_step(step);
                auto node_length = graph.get_length(handle);
                auto node_end = offset + node_length;
                // remove non-active ranges and call callbacks
                active_ranges.erase(
                    std::remove_if(active_ranges.begin(),
                                   active_ranges.end(),
                                   [&](auto& r) {
                                       if (r.first->end.offset <= offset) {
                                           func(*r.first,
                                                (double)r.second /
                                                (double)(r.first->end.offset
                                                         - r.first->begin.offset));
                                           return true;
                                       } else {
                                           return false;
                                       }
                                   }),
                    active_ranges.end());
                // find ranges that start in this node, and add them to active ranges
                while (range_itr != ranges_end
                       && (*range_itr)->begin.offset < node_end
                       && (*range_itr)->begin.offset >= offset) {
                    active_ranges.push_back(std::make_pair(*range_itr++, 0));
                }
                // compute coverage on node
                auto& d = depths[graph.get_id(handle) - shift];
                // add coverage to active ranges
                for (auto& r : active_ranges) {
                    auto l = node_length;
                    if (r.first->begin.offset > offset) {
                        l -= r.first->begin.offset - offset;
                    }
                    if (r.first->end.offset < node_end) {
                        l -= node_end - r.first->end.offset;
                    }
                    r.second += (d * l);
                }
                // to node end
                offset = node_end;
            });
        // clear out any ranges that are still active at the end
        for (auto& r : active_ranges) {
            func(*r.first,
                 (double)r.second /
                 (double)(r.first->end.offset
                          - r.first->begin.offset));
        }
    }
}

}
}
