#include "bin_path_info.hpp"

namespace odgi {
namespace algorithms {

void bin_path_coverage(const PathHandleGraph& graph,
                       const std::string& prefix_delimiter,
                       const uint64_t& num_bins,
                       std::vector<std::pair<std::string, std::vector<double>>>& table) {
    // the graph must be compacted for this to work
    std::vector<uint64_t> position_map(graph.get_node_count()+1);
    std::vector<std::pair<uint64_t, uint64_t>> contacts;
    uint64_t len = 0;
    graph.for_each_handle([&](const handle_t& h) {
            position_map[number_bool_packing::unpack_number(h)] = len;
            uint64_t hl = graph.get_length(h);
            len += hl;
        });
    position_map[position_map.size()-1] = len;
    // find the prefix order based on existing set
    std::unordered_set<std::string> seen_prefixes;
    std::unordered_map<std::string, uint64_t> path_prefix_rank;
    auto get_path_prefix = [&](const path_handle_t& p) -> std::string {
        std::string path_name = graph.get_path_name(p);
        if (prefix_delimiter.empty()) {
            return path_name;
        } else {
            return path_name.substr(0, path_name.find(prefix_delimiter));
        }
    };
    graph.for_each_path_handle([&](const path_handle_t& p) {
            std::string path_prefix = get_path_prefix(p);
            if (!seen_prefixes.count(path_prefix)) {
                path_prefix_rank[path_prefix] = seen_prefixes.size();
                seen_prefixes.insert(path_prefix);
            }
        });
    std::vector<std::string> path_prefix_order(seen_prefixes.size());
    for (auto& p : path_prefix_rank) {
        path_prefix_order[p.second] = p.first;
    }
    // resize the aggregation matrix
    table.resize(path_prefix_rank.size());
    auto path_prefix_order_itr = path_prefix_order.begin();
    for (auto& row : table) {
        row.first = *path_prefix_order_itr++;
        row.second.resize(num_bins);
    }
    graph.for_each_path_handle([&](const path_handle_t& path) {
            auto& path_rank = path_prefix_rank[get_path_prefix(path)];
            auto& row = table[path_rank].second;
            // walk the path and aggregate
            graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
                    handle_t h = graph.get_handle_of_step(occ);
                    uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                    uint64_t hl = graph.get_length(h);
                    // make contects for the bases in the node
                    //std::cerr << p << std::endl;
                    for (uint64_t k = 0; k < hl; ++k) {
                        double x = std::floor((double)(p+k)/(double)len * num_bins);
                        ++row[x];
                    }
                });
        });
    double bp_per_bin = (double)len / (double)num_bins;
    for (auto& row : table) {
        for (auto& v : row.second) {
            v /= bp_per_bin;
        }
    }
}

}
}
