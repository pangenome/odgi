#include "bin_path_info.hpp"

namespace odgi {
namespace algorithms {

void bin_path_info(const PathHandleGraph& graph,
                   const std::string& prefix_delimiter,
                   const std::function<void(const std::string&, const hash_map<uint64_t, algorithms::path_info_t>&)>& handle_path,
                   uint64_t num_bins,
                   uint64_t bin_width) {
    // the graph must be compacted for this to work
    std::vector<uint64_t> position_map(graph.get_node_count()+1);
    std::vector<std::pair<uint64_t, uint64_t>> contacts;
    uint64_t len = 0;
    graph.for_each_handle([&](const handle_t& h) {
            position_map[number_bool_packing::unpack_number(h)] = len;
            uint64_t hl = graph.get_length(h);
            len += hl;
        });
    if (!num_bins) {
        num_bins = len / bin_width + (len % bin_width ? 1 : 0);
    } else if (!bin_width) {
        bin_width = len / num_bins;
        num_bins = len / bin_width + (len % bin_width ? 1 : 0);
    }
    position_map[position_map.size()-1] = len;
    hash_map<path_handle_t, uint64_t> path_length;
    graph.for_each_path_handle([&](const path_handle_t& path) {
            hash_map<uint64_t, path_info_t> table;
            // walk the path and aggregate
            uint64_t path_pos = 0;
            int64_t last_bin = -1; // flag meaning "null bin"
            uint64_t last_pos_in_bin = 0;
            bool last_is_rev = false;
            graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
                    handle_t h = graph.get_handle_of_step(occ);
                    bool is_rev = graph.get_is_reverse(h);
                    uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                    uint64_t hl = graph.get_length(h);
                    // detect bin crossings
                    // make contects for the bases in the node
                    for (uint64_t k = 0; k < hl; ++k) {
                        int64_t curr_bin = (p+k) / bin_width;
                        uint64_t curr_pos_in_bin = (p+k) - (curr_bin * bin_width);
                        if (curr_bin != last_bin) {
                            // bin cross!
                            table[curr_bin].begins.push_back({last_bin,last_pos_in_bin,curr_pos_in_bin,path_pos,is_rev});
                            if (last_bin >= 0) {
                                table[last_bin].ends.push_back({curr_bin,curr_pos_in_bin,last_pos_in_bin,path_pos-1,last_is_rev});
                            }
                        }
                        ++table[curr_bin].mean_cov;
                        if (is_rev) {
                            ++table[curr_bin].mean_inv;
                        }
                        table[curr_bin].mean_pos += path_pos++;
                        last_bin = curr_bin;
                        last_is_rev = is_rev;
                        last_pos_in_bin = curr_pos_in_bin;
                    }
                });
            // record the last one
            table[last_bin].ends.push_back({-1,0,last_pos_in_bin,path_pos-1,last_is_rev});
            uint64_t path_length = path_pos;
            for (auto& entry : table) {
                auto& v = entry.second;
                v.mean_inv /= (v.mean_cov ? v.mean_cov : 1);
                v.mean_cov /= bin_width;
                v.mean_pos /= bin_width * path_length * v.mean_cov;
            }
            handle_path(graph.get_path_name(path), table);
        });
}

}
}
