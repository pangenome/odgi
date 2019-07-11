#include "id_ordered_paths.hpp"

namespace odgi {
namespace algorithms {

std::vector<path_handle_t> id_ordered_paths(const PathHandleGraph& g, bool avg, bool rev) {
    return prefix_and_id_ordered_paths(g, "", avg, rev);
}

std::vector<path_handle_t> prefix_and_id_ordered_paths(const PathHandleGraph& g, const std::string& prefix_delimiter, bool avg, bool rev) {
    // find the prefix order based on existing set
    std::unordered_set<std::string> seen_prefixes;
    std::unordered_map<std::string, uint64_t> initial_path_prefix_rank;
    auto get_path_prefix = [&](const path_handle_t& p) -> std::string {
        if (prefix_delimiter.empty()) {
            return "";
        } else {
            std::string path_name = g.get_path_name(p);
            return path_name.substr(0, path_name.find(prefix_delimiter));
        }
    };
    g.for_each_path_handle([&](const path_handle_t& p) {
            std::string path_prefix = get_path_prefix(p);
            if (!seen_prefixes.count(path_prefix)) {
                initial_path_prefix_rank[path_prefix] = seen_prefixes.size();
                seen_prefixes.insert(path_prefix);
            }
        });
    std::vector<std::string> path_prefix_order(seen_prefixes.size());
    for (auto& p : initial_path_prefix_rank) {
        path_prefix_order[p.second] = p.first;
    }
    // accumulate the paths by their minimum id per prefix
    std::unordered_map<std::string, std::vector<std::pair<double, path_handle_t> > > paths;
    g.for_each_path_handle([&](const path_handle_t& p) {
            // get the min id
            if (avg) {
                double sum_id = 0;
                uint64_t step_count = 0;
                g.for_each_step_in_path(p, [&](const step_handle_t& occ) {
                        sum_id += (double)g.get_id(g.get_handle_of_step(occ));
                        ++step_count;
                    });
                std::string path_prefix = get_path_prefix(p);
                paths[path_prefix].push_back(make_pair(sum_id/(double)step_count, p));
            } else {
                double min_id = std::numeric_limits<double>::max();
                g.for_each_step_in_path(p, [&](const step_handle_t& occ) {
                        min_id = std::min((double)g.get_id(g.get_handle_of_step(occ)), min_id);
                    });
                std::string path_prefix = get_path_prefix(p);
                paths[path_prefix].push_back(make_pair(min_id, p));
            }
        });
    for (auto& b : paths) {
        auto& paths_by = b.second;
        std::sort(paths_by.begin(), paths_by.end(), [&](const std::pair<double, path_handle_t>& a,
                                                        const std::pair<double, path_handle_t>& b) {
                      return a.first < b.first || a.first == b.first && as_integer(a.second) < as_integer(b.second);
                  });
        if (rev) {
            std::reverse(paths_by.begin(), paths_by.end());
        }
    }
    std::vector<path_handle_t> order; order.reserve(g.get_path_count());
    for (auto& path_prefix : path_prefix_order) {
        for (auto& p : paths[path_prefix]) {
            order.push_back(p.second);
        }
    }
    return order;
}

}
}
