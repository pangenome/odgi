#include "bin_path_info.hpp"

// #define  debug_bin_path_info

namespace odgi {
    namespace algorithms {

        template<int first, int second>
        void updatePair(std::pair<uint64_t, uint64_t> &p, uint64_t nucleotide_count) {
#ifdef debug_bin_path_info
            std::cerr << "UPDATE_PAIR: " << "<" << std::get<first>(p) << "," << std::get<second>(p) << ">" << std::endl;
            std::cerr << "NUC: " << nucleotide_count << std::endl;
#endif
            if (std::get<second>(p) == 0) {
                std::get<second>(p) = std::get<first>(p);

            }
            std::get<first>(p) = nucleotide_count;
#ifdef debug_bin_path_info
            std::cerr << "AFTER_UPDATE_PAIR: " << "<" << std::get<first>(p) << "," << std::get<second>(p) << ">"
                      << std::endl;
#endif
        }

        void bin_path_info(const PathHandleGraph &graph,
                           const std::string &prefix_delimiter,
                           const std::function<void(const uint64_t &, const uint64_t &)> &handle_header,
                           const std::function<void(const std::string &,
                                                    const std::vector<std::pair<uint64_t, uint64_t>> &,
                                                    const std::map<uint64_t, algorithms::path_info_t> &)> &handle_path,
                           const std::function<void(const uint64_t &, const std::string &)> &handle_sequence,
                           uint64_t num_bins,
                           uint64_t bin_width,
                           bool drop_gap_links,
                           bool progress) {
            // the graph must be compacted for this to work
            std::vector<uint64_t> position_map(graph.get_node_count() + 1);
            uint64_t len = 0;
            std::string graph_seq;
            graph.for_each_handle([&](const handle_t &h) {
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
            position_map[position_map.size() - 1] = len;
            // write header
            handle_header(len, bin_width);
            // collect bin sequences
            for (uint64_t i = 0; i < num_bins; ++i) {
                handle_sequence(i + 1, graph_seq.substr(i * bin_width, bin_width));
            }
            graph_seq.clear(); // clean up
            std::unordered_map<path_handle_t, uint64_t> path_length;
            uint64_t gap_links_removed = 0;
            uint64_t total_links = 0;
            std::unique_ptr<progress_meter::ProgressMeter> progress_meter;
            if (progress) {
                progress_meter = std::make_unique<progress_meter::ProgressMeter>(
                        graph.get_path_count(), "[odgi::bin] bin_path_info:");
            }
            graph.for_each_path_handle([&](const path_handle_t &path) {
                std::vector<std::pair<uint64_t, uint64_t>> links;
                std::map<uint64_t, path_info_t> bins;
                // walk the path and aggregate
                uint64_t path_pos = 0;
                int64_t last_bin = 0; // flag meaning "null bin"
                uint64_t last_pos_in_bin = 0;
                uint64_t nucleotide_count = 0;
                bool last_is_rev = false;
                graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                    handle_t h = graph.get_handle_of_step(occ);
                    bool is_rev = graph.get_is_reverse(h);
                    uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                    uint64_t hl = graph.get_length(h);
                    // detect bin crossings
                    // make contects for the bases in the node
                    for (uint64_t k = 0; k < hl; ++k) {
                        int64_t curr_bin = (p + k) / bin_width + 1;
                        uint64_t curr_pos_in_bin = (p + k) - (curr_bin * bin_width);
                        if (curr_bin != last_bin && std::abs(curr_bin - last_bin) > 1 || last_bin == 0) {
                            // bin cross!
                            links.push_back(std::make_pair(last_bin, curr_bin));
                        }
                        ++bins[curr_bin].mean_depth;
                        if (is_rev) {
                            ++bins[curr_bin].mean_inv;
                        }
                        bins[curr_bin].mean_pos += path_pos++;
                        nucleotide_count += 1;
                        if ((bins[curr_bin].ranges.size() == 0) ||
                            ((nucleotide_count - bins[curr_bin].ranges.back().second) > 1 &&
                             (nucleotide_count - bins[curr_bin].ranges.back().first) > 1) ||
                            (is_rev != last_is_rev)) {
                            std::pair<uint64_t, uint64_t> p = std::make_pair(0, 0);
                            if (is_rev) {
                                std::get<0>(p) = nucleotide_count;
                            } else {
                                std::get<1>(p) = nucleotide_count;
                            }
                            bins[curr_bin].ranges.push_back(p);
#ifdef debug_bin_path_info
                            std::cerr << "PUSHED PAIR: " << "<" << std::get<0>(p) << "," << std::get<1>(p) << ">"
                                      << std::endl;
#endif
                        } else {
                            std::pair<uint64_t, uint64_t> &p = bins[curr_bin].ranges.back();
                            if (is_rev) {
                                updatePair<0, 1>(p, nucleotide_count);
                            }
                            else {
                                updatePair<1, 0>(p, nucleotide_count);
                            }
                        }
                        last_bin = curr_bin;
                        last_is_rev = is_rev;
                        last_pos_in_bin = curr_pos_in_bin;
                    }
                });
                links.push_back(std::make_pair(last_bin, 0));
                uint64_t path_length = path_pos;
                uint64_t end_nucleotide = nucleotide_count;
                for (auto &entry : bins) {
                    auto &v = entry.second;
                    v.mean_inv /= (v.mean_depth ? v.mean_depth : 1);
                    v.mean_depth /= bin_width;
                    v.mean_pos /= bin_width * path_length * v.mean_depth;
                }

                if (drop_gap_links) {
                    std::vector<uint64_t> bin_ids;
                    for (const auto &entry: bins) {
                        bin_ids.push_back(entry.first);
                    }
                    std::sort(bin_ids.begin(), bin_ids.end());
                    total_links += links.size();

                    uint64_t fill_pos = 0;

                    for (uint64_t i = 0; i < links.size(); ++i) {
                        auto link = links[i];

                        if (link.first == 0 || link.second == 0)
                            continue;

                        if (link.first > link.second) {
                            links[fill_pos++] = link;
                            continue;
                        }

                        auto left_it = std::lower_bound(bin_ids.begin(), bin_ids.end(), link.first + 1);
                        auto right_it = std::lower_bound(bin_ids.begin(), bin_ids.end(), link.second);
                        if (right_it > left_it) {
                            links[fill_pos++] = link;
                        }
                    }

                    gap_links_removed += links.size() - fill_pos;
                    links.resize(fill_pos);
                }

                handle_path(graph.get_path_name(path), links, bins);

                if (progress) {
                    progress_meter->increment(1);
                }
            });

            if (progress) {
                progress_meter->finish();
            }

            if (drop_gap_links) {
                uint64_t path_count = graph.get_path_count();

                std::cerr << std::setprecision(4) << "[odgi::bin_path_info] Gap links removed: " << (100.0 *  ((double)gap_links_removed / (double)total_links))
                          << "%, that is " << gap_links_removed << " gap links (" << path_count << " path start links + "
                          << path_count << " path end links + " << (gap_links_removed - path_count * 2) << " inner gap links) of "
                          << total_links << " total links" << std::endl;
            }
        }

    }
}
