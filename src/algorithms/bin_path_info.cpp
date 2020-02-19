#include <stdint.h>
#include <cstdint>
#include <random>
#include <chrono>
#include <string>
#include <fstream>

#include "bin_path_info.hpp"

namespace odgi {
namespace algorithms {

    template<typename T> size_t get_nbytes(const T &c) {
        std::ofstream ofs("/dev/null");
        return c.serialize(ofs);
    }

void bin_path_info(const PathHandleGraph& graph,
                   const std::string& prefix_delimiter,
                   const std::function<void(const std::string&,
                                            const std::vector<std::pair<uint64_t, uint64_t>>&,
                                            const std::map<uint64_t, algorithms::path_info_t>&)>& handle_path,
                   const std::function<void(const uint64_t&, const std::string&)>& handle_sequence,
                   // const std::function<void(const std::map<std::string&, std::vector<uint64_t&>>) handle_index, // TODO @ekg I did not read up yet about pointers again. Would "&" be appropriate here?
                   uint64_t num_bins,
                   uint64_t bin_width,
                   bool index_json) {
    // the graph must be compacted for this to work
    std::vector<uint64_t> position_map(graph.get_node_count()+1);
    uint64_t len = 0;
    std::string graph_seq;
    graph.for_each_handle([&](const handle_t& h) {
            position_map[number_bool_packing::unpack_number(h)] = len;
            uint64_t hl = graph.get_length(h);
            graph_seq.append(graph.get_sequence(h));
            len += hl;
        });
    if (!num_bins) {
        num_bins = len / bin_width + (len % bin_width ? 1 : 0);
    } else if (!bin_width) {
        bin_width = len / num_bins;
        num_bins = len / bin_width + (len % bin_width ? 1 : 0);
    }
    position_map[position_map.size()-1] = len;
    // collect bin sequences
    for (uint64_t i = 0; i < num_bins; ++i) {
        handle_sequence(i+1, graph_seq.substr(i*bin_width, bin_width));
    }
    graph_seq.clear(); // clean up
    std::unordered_map<path_handle_t, uint64_t> path_length;

    // path position index for the bins
    std::map<std::string, std::vector<uint64_t>> bin_index;

    graph.for_each_path_handle([&](const path_handle_t& path) {
            std::vector<std::pair<uint64_t, uint64_t>> links;
            std::map<uint64_t, path_info_t> bins;

            // walk the path and aggregate
            uint64_t path_pos = 0;
            int64_t last_bin = 0; // flag meaning "null bin"
            uint64_t last_pos_in_bin = 0;
            bool last_is_rev = false;

            // this is a dynamic integer vector
            std::vector<uint64_t> path_bin_vec;

            graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {

                    handle_t h = graph.get_handle_of_step(occ);
                    bool is_rev = graph.get_is_reverse(h);
                    uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                    uint64_t hl = graph.get_length(h); // TODO @ekg Is this the length of the current handle of the current path in nucleotides? Or the number of nodes?
                    // detect bin crossings
                    // make contects for the bases in the node
                    for (uint64_t k = 0; k < hl; ++k) {
                        int64_t curr_bin = (p+k) / bin_width + 1;
                        uint64_t curr_pos_in_bin = (p+k) - (curr_bin * bin_width);
                        if (curr_bin != last_bin && std::abs(curr_bin-last_bin) > 1 || last_bin == 0) {
                            // bin cross!
                            links.push_back(std::make_pair(last_bin,curr_bin));
                        }
                        ++bins[curr_bin].mean_cov;
                        if (is_rev) {
                            ++bins[curr_bin].mean_inv;
                        }

                        path_pos++;
                        // FIXME at the current path position we add the current bin we are in into the compressed integer vector
                        path_bin_vec.push_back(curr_bin);

                        bins[curr_bin].mean_pos += path_pos;
                        last_bin = curr_bin;
                        last_is_rev = is_rev;
                        last_pos_in_bin = curr_pos_in_bin;
                    }
                });
            links.push_back(std::make_pair(last_bin,0));
            uint64_t path_length = path_pos;
            for (auto& entry : bins) {
                auto& v = entry.second;
                v.mean_inv /= (v.mean_cov ? v.mean_cov : 1);
                v.mean_cov /= bin_width;
                v.mean_pos /= bin_width * path_length * v.mean_cov;
            }
            std::string path_name = graph.get_path_name(path);

            std::cout << "compression_method" << '\t' << "path_name" << '\t' << "num_bytes" << "\n";
            // TODO Get the template from  https://github.com/dnbaker/sdsl-vec-test/blob/master/src/sdsl-vec-test.cpp so that I can reuse code.

            // create a compressed integer vector via sdsl-lite
            sdsl::enc_vector<sdsl::coder::elias_gamma> elias_gamma_path_bin_vec(path_bin_vec);
            size_t nb = get_nbytes(elias_gamma_path_bin_vec);
            std::cout << "elias_gamma" << '\t' << path_name << '\t' << nb  << "\n";
            // nb = get_nbytes(path_bin_vec); FIXME @ekg This does not work, any ideas? In https://github.com/dnbaker/sdsl-vec-test/blob/master/src/sdsl-vec-test.cpp it seems to work that way.
            // std::cout << "uncompressed" << '\t' << path_name << '\t' << nb  << "\n"; TODO
            sdsl::enc_vector<sdsl::coder::elias_delta> elias_delta_path_bin_vec(path_bin_vec);
            nb = get_nbytes(elias_delta_path_bin_vec);
            std::cout << "elias_delta" << '\t' << path_name << '\t' << nb  << "\n";

            sdsl::enc_vector<sdsl::coder::fibonacci> fibonacci_path_bin_vec(path_bin_vec);
            nb = get_nbytes(fibonacci_path_bin_vec);
            std::cout << "fibonacci" << '\t' << path_name << '\t' << nb  << "\n";

            sdsl::vlc_vector<sdsl::coder::elias_gamma> vlc_elias_gamma_path_bin_vec(path_bin_vec);
            nb = get_nbytes(vlc_elias_gamma_path_bin_vec);
            std::cout << "vlc_elias_gamma" << '\t' << path_name << '\t' << nb  << "\n";

            // FIXME @ekg I want to add the compressed_integer_vector of the current path to the map. However, the variable won't be reachable any more. Any workarounds?
            // bin_index[path_name] = compressed_integer_vector
            // if (index_json) {
            //     handle_index(bin_index); // TODO Write the handle_index function
            // }
            handle_path(graph.get_path_name(path), links, bins);
        });
}

}
}
