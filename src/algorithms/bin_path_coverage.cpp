#include "bin_path_coverage.hpp"

// #define  debug_bin_path_coverage

namespace odgi {
    namespace algorithms {

        void bin_path_coverage(const PathHandleGraph &graph,
                           uint64_t num_bins,
                           uint64_t bin_width,
                           bool progress) {
            // the graph must be compacted for this to work
            std::vector<uint64_t> position_map(graph.get_node_count() + 1);
            uint64_t graph_len = 0;
            graph.for_each_handle([&](const handle_t &h) {
                position_map[number_bool_packing::unpack_number(h)] = graph_len;
                uint64_t hl = graph.get_length(h);
                graph_len += hl;
            });
            if (!num_bins) {
                num_bins = graph_len / bin_width + (graph_len % bin_width ? 1 : 0);
            } else if (!bin_width) {
                bin_width = graph_len / num_bins;
                num_bins = graph_len / bin_width + (graph_len % bin_width ? 1 : 0);
            }
            position_map[position_map.size() - 1] = graph_len;
            std::cerr <<"[odgi::bin] bin_path_coverage: pangenome length: " << graph_len << std::endl;
            std::vector<uint64_t> paths_per_bin(num_bins);
            // first pass: collect #nucleotides, fill in_all_bins_bv, unique_bins_bv, unique_bins_touched_bv
            graph.for_each_path_handle([&](const path_handle_t &path) {
                std::vector<bool> paths_per_bin_bv(num_bins);
                graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                    handle_t h = graph.get_handle_of_step(occ);
                    uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                    uint64_t hl = graph.get_length(h);
                    for (uint64_t k = 0; k < hl; ++k) {
                        int64_t curr_bin = (p + k) / bin_width + 1;
                        uint64_t curr_pos_in_bin = (p + k) - (curr_bin * bin_width);
                        if (!paths_per_bin_bv[curr_bin - 1]) {
                            paths_per_bin[curr_bin - 1] += 1;
                            paths_per_bin_bv[curr_bin - 1] = true;
                        }
                    }
                });
            });
#ifdef debug_bin_path_coverage
            for (int i = 0; i < paths_per_bin.size(); i++) {
                std::cerr << paths_per_bin[i] << std::endl;
            }
#endif

            // write header of table to stdout, collect final bins
            std::vector<uint64_t> final_bin_ids = std::vector<uint64_t>();
            const uint64_t graph_paths = graph.get_path_count();
            std::cout << "path_name";
            for (uint64_t i = 0; i < num_bins; i++) {
                uint64_t paths_in_bin = paths_per_bin[i];
                // filter out bins that are unique (== 1) or that are present in all paths (== graph_paths)
                // 1-based bin identifiers
                if (!(paths_in_bin == 1 || paths_in_bin == graph_paths)) {
                    std::cout << "\t" << (i + 1);
                    final_bin_ids.push_back(i + 1);
                }
            }
            std::cout << std::endl;

            // for all paths, for each node and nucleotide in path --> better unordered_map https://stackoverflow.com/questions/1939953/how-to-find-if-a-given-key-exists-in-a-c-stdmap
            // map<uint64_t, double>: ++ when we cover the bin
            // as it is done so far
            std::unique_ptr<progress_meter::ProgressMeter> progress_meter;
            if (progress) {
                progress_meter = std::make_unique<progress_meter::ProgressMeter>(
                        graph.get_path_count(), "[odgi::bin] bin_path_coverage:");
            }

            graph.for_each_path_handle([&](const path_handle_t &path) {
                std::unordered_map<uint64_t, double> bin_coverages;
                graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                    handle_t h = graph.get_handle_of_step(occ);
                    bool is_rev = graph.get_is_reverse(h);
                    uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                    uint64_t hl = graph.get_length(h);
                    for (uint64_t k = 0; k < hl; ++k) {
                        int64_t curr_bin = (p + k) / bin_width + 1;
                        uint64_t curr_pos_in_bin = (p + k) - (curr_bin * bin_width);
                        // update mean coverage
                        bin_coverages[curr_bin] += 1;
                    }
                });
                for (auto &entry : bin_coverages) {
                    auto &v = entry.second;
                    auto &k = entry.first;
                    bin_coverages[k] = bin_coverages[k] / bin_width;
                }

                // write to std::cout, only if we have !in_all_bins_bv[idx] && !unique_bins_bv[idx] && auto it = m.find("f"); if (it != m.end()) {/*Use it->second*/}.
                std::cout << graph.get_path_name(path);
                // this could be optimized by only iterating over the bins that we want to print, this would require some pre-calculations, might be memory expensive
                for (auto &bin_id : final_bin_ids) {
                    auto it = bin_coverages.find(bin_id);
                    if (it != bin_coverages.end()) {
                        uint64_t cov = it->second;
                        // cap coverage by 255
                        if (cov > 255) {
                            cov = 255;
                        }
                        std::cout << "\t" << cov;
                    } else {
                        std::cout << "\t" << 0;
                    }
                }
                std::cout << std::endl;
                // write header
                if (progress) {
                    progress_meter->increment(1);
                }
            });

            if (progress) {
                progress_meter->finish();
            }
        }
    }
}
