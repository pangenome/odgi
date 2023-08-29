#include "extract.hpp"

namespace odgi {
    namespace algorithms {

        void add_full_paths_to_component(const graph_t &source, graph_t &component, const uint64_t num_threads) {
            // Search paths in parallel
            atomicbitvector::atomic_bv_t take_source_path(source.get_path_count());

            component.for_each_handle([&](const handle_t &h) {
                const handlegraph::nid_t id = component.get_id(h);

                if (source.has_node(id)) {
                    handle_t source_handle = source.get_handle(id);

                    source.for_each_step_on_handle(source_handle, [&](const step_handle_t &source_step) {
                        const uint64_t source_path_rank = as_integer(source.get_path_handle_of_step(source_step)) - 1;
                        take_source_path.set(source_path_rank);
                    });
                }
            }, true);

            // Create paths
            std::vector<path_handle_t> taken_source_paths;
            source.for_each_path_handle([&](const path_handle_t source_path) {
                const uint64_t source_path_rank = as_integer(source_path) - 1;

                if (take_source_path.test(source_path_rank)) {
                    component.create_path_handle(source.get_path_name(source_path), source.get_is_circular(source_path));
                    taken_source_paths.push_back(source_path);
                }
            });

            // Fill paths in parallel
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
            for (auto& source_path : taken_source_paths) {
                const path_handle_t path_handle = component.get_path_handle(source.get_path_name(source_path));

                for (handle_t handle : source.scan_path(source_path)) {
                    component.append_step(path_handle, component.get_handle(source.get_id(handle),
                                                                                source.get_is_reverse(handle)));
                }
            }
        }

        // Create a subpath name
        std::string make_path_name(const string &path_name, size_t offset, size_t end_offset) {
            std::string out_name = path_name + ":" + std::to_string(offset);
            if (end_offset > offset) {
                out_name += "-" + std::to_string(end_offset);
            }
            return out_name;
        }

        path_handle_t create_subpath(graph_t &subgraph, const string &subpath_name, const bool is_circular) {
            // The function assumes that every path is new and unique
            // if (subgraph.has_path(subpath_name)) {
            //     subgraph.destroy_path(subgraph.get_path_handle(subpath_name));
            // }
            return subgraph.create_path_handle(subpath_name, is_circular);
        };

        void add_subpaths_to_subgraph(const graph_t &source, const std::vector<path_handle_t> source_paths,
                                      graph_t &subgraph, const uint64_t num_threads,
                                      const std::string &progress_message) {
            const bool show_progress = !progress_message.empty();

            std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress;
            if (show_progress) {
                progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                        source_paths.size() * 3, progress_message);
            }

            std::vector<std::vector<std::pair<uint64_t, uint64_t>>> subpath_ranges;
            subpath_ranges.resize(source_paths.size());

            // Search subpaths in parallel
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
            for (uint64_t path_rank = 0; path_rank < source_paths.size(); ++path_rank) {
                auto &source_path_handle = source_paths[path_rank];
                const std::string path_name = source.get_path_name(source_path_handle);

                uint64_t walked = 0;
                bool first_node = true;

                uint64_t start, end;
                source.for_each_step_in_path(source_path_handle, [&](const step_handle_t &step) {
                    handle_t source_handle = source.get_handle_of_step(step);
                    uint64_t source_length = source.get_length(source_handle);
                    handlegraph::nid_t source_id = source.get_id(source_handle);

                    if (subgraph.has_node(source_id)) {
                        if (first_node) {
                            first_node = false;
                            start = walked;
                        }

                        end = walked + source_length;
                    } else if (!first_node) {
                        subpath_ranges[path_rank].push_back({start, end});
                        first_node = true;
                    }

                    walked += source_length;
                });

                // last subpath
                if (!first_node) {
                    subpath_ranges[path_rank].push_back({start, end});
                }

                if (show_progress) {
                    progress->increment(1);
                }
            }

            // Create subpaths
            for (uint64_t path_rank = 0; path_rank < source_paths.size(); ++path_rank) {
                const auto &source_path_handle = source_paths[path_rank];
                const std::string path_name = source.get_path_name(source_path_handle);

                for (auto subpath_range : subpath_ranges[path_rank]) {
                    create_subpath(subgraph, make_path_name(path_name, subpath_range.first, subpath_range.second),
                                   source.get_is_circular(source_path_handle));
                }

                if (show_progress) {
                    progress->increment(1);
                }
            }

            // Fill subpaths in parallel
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
            for (uint64_t path_rank = 0; path_rank < source_paths.size(); ++path_rank) {
                if (!subpath_ranges[path_rank].empty()) {
                    const auto &source_path_handle = source_paths[path_rank];
                    const std::string path_name = source.get_path_name(source_path_handle);

                    // The path ranges are sorted by coordinates by design
                    uint64_t range_rank = 0;
                    path_handle_t subpath_handle = subgraph.get_path_handle(
                            make_path_name(path_name, subpath_ranges[path_rank][0].first,
                                           subpath_ranges[path_rank][0].second)
                    );

                    uint64_t walked = 0;
                    source.for_each_step_in_path(source_path_handle, [&](const step_handle_t &step) {
                        if (range_rank < subpath_ranges[path_rank].size()) {
                            const handle_t source_handle = source.get_handle_of_step(step);

                            if (walked >= subpath_ranges[path_rank][range_rank].first &&
                                walked <= subpath_ranges[path_rank][range_rank].second) {
                                subgraph.append_step(
                                        subpath_handle,
                                        subgraph.get_handle(source.get_id(source_handle),
                                                            source.get_is_reverse(source_handle))
                                );
                            }

                            walked += source.get_length(source_handle);
                            if (walked >= subpath_ranges[path_rank][range_rank].second) {
                                ++range_rank;
                                if (range_rank < subpath_ranges[path_rank].size()) {
                                    subpath_handle = subgraph.get_path_handle(
                                            make_path_name(path_name, subpath_ranges[path_rank][range_rank].first,
                                                           subpath_ranges[path_rank][range_rank].second)
                                    );
                                }//else not other subpath ranges for this path
                            }
                        }
                    });
                }

                if (show_progress) {
                    progress->increment(1);
                }
            }

            if (show_progress) {
                progress->finish();
            }
        }

        void extract_path_range(const graph_t &source, path_handle_t path_handle, int64_t start, int64_t end,
                                graph_t &subgraph) {
            algorithms::for_handle_in_path_range(
                    source, path_handle, start, end,
                    [&](const handle_t& cur_handle) {
                        const nid_t id = source.get_id(cur_handle);
                        if (!subgraph.has_node(id)) {
                            subgraph.create_handle(
                                    source.get_sequence(
                                            source.get_is_reverse(cur_handle) ? source.flip(cur_handle) : cur_handle),
                                            id);
                        }
                    });
        }

        void for_handle_in_path_range(const graph_t &source, path_handle_t path_handle, int64_t start, int64_t end,
                                      const std::function<void(const handle_t&)>& lambda) {
            uint64_t walked = 0;
            const auto path_end = source.path_end(path_handle);
            for (step_handle_t cur_step = source.path_begin(path_handle);
                 cur_step != path_end && walked < end; cur_step = source.get_next_step(cur_step)) {
                const handle_t cur_handle = source.get_handle_of_step(cur_step);
                walked += source.get_length(cur_handle);
                if (walked > start) {
                    lambda(cur_handle);
                }
            }
        }

        /// We can accumulate a subgraph without accumulating all the edges between its nodes
        /// this helper ensures that we get the full set
        void add_connecting_edges_to_subgraph(const graph_t &source, graph_t &subgraph,
                                              const std::string &progress_message) {
            const bool show_progress = !progress_message.empty();

            std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress;
            if (show_progress) {
                progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                        subgraph.get_node_count(), progress_message);
            }

            subgraph.for_each_handle([&](const handle_t &handle) {
                const nid_t id = subgraph.get_id(handle);
                const handle_t source_handle = source.get_handle(id, subgraph.get_is_reverse(handle));
                source.follow_edges(source_handle, false, [&](const handle_t &next) {
                    const nid_t next_id = source.get_id(next);
                    if (subgraph.has_node(next_id)) {
                        handle_t subgraph_next = subgraph.get_handle(next_id, source.get_is_reverse(next));
                        if (!subgraph.has_edge(handle, subgraph_next)) {
                            subgraph.create_edge(handle, subgraph_next);
                        }
                    }
                });
                source.follow_edges(source_handle, true, [&](const handle_t &prev) {
                    const nid_t prev_id = source.get_id(prev);
                    if (subgraph.has_node(prev_id)) {
                        const handle_t subgraph_prev = subgraph.get_handle(prev_id, source.get_is_reverse(prev));
                        if (!subgraph.has_edge(subgraph_prev, handle)) {
                            subgraph.create_edge(subgraph_prev, handle);
                        }
                    }
                });

                if (show_progress) {
                    progress->increment(1);
                }
            }, false);

            if (show_progress) {
                progress->finish();
            }
        }

        void
        expand_subgraph_by_steps(const graph_t &source, graph_t &subgraph, const uint64_t &steps, bool forward_only) {
            std::vector<handle_t> curr_handles;
            subgraph.for_each_handle([&](const handle_t &h) {
                curr_handles.push_back(h);
            });
            for (uint64_t i = 0; i < steps && !curr_handles.empty(); ++i) {
                std::vector<handle_t> next_handles;
                for (auto &h : curr_handles) {
                    const handle_t old_h = source.get_handle(subgraph.get_id(h));
                    source.follow_edges(old_h, false, [&](const handle_t &c) {
                        if (!subgraph.has_node(source.get_id(c))) {
                            const handle_t x = subgraph.create_handle(
                                    source.get_sequence(source.get_is_reverse(c) ? source.flip(c) : c),
                                    source.get_id(c));
                            next_handles.push_back(x);
                        }
                    });
                    if (!forward_only) {
                        source.follow_edges(old_h, true, [&](const handle_t &c) {
                            if (!subgraph.has_node(source.get_id(c))) {
                                const handle_t x = subgraph.create_handle(
                                        source.get_sequence(source.get_is_reverse(c) ? source.flip(c) : c),
                                        source.get_id(c));
                                next_handles.push_back(x);
                            }
                        });
                    }
                }
                curr_handles = std::move(next_handles);
            }
        }

        void
        expand_subgraph_by_length(const graph_t &source, graph_t &subgraph, const uint64_t &length, bool forward_only) {
            uint64_t accumulated_length = 0;
            std::vector<handle_t> curr_handles;
            subgraph.for_each_handle([&](const handle_t &h) {
                curr_handles.push_back(h);
            });
            while (accumulated_length < length && !curr_handles.empty()) {
                std::vector<handle_t> next_handles;
                for (auto &h : curr_handles) {
                    const handle_t old_h = source.get_handle(subgraph.get_id(h));
                    source.follow_edges(old_h, false, [&](const handle_t &c) {
                        if (!subgraph.has_node(source.get_id(c))) {
                            const handle_t x = subgraph.create_handle(
                                    source.get_sequence(source.get_is_reverse(c) ? source.flip(c) : c),
                                    source.get_id(c));
                            next_handles.push_back(x);
                            accumulated_length += subgraph.get_length(x);
                        }
                    });
                    if (!forward_only) {
                        source.follow_edges(old_h, true, [&](const handle_t &c) {
                            if (!subgraph.has_node(source.get_id(c))) {
                                const handle_t x = subgraph.create_handle(
                                        source.get_sequence(source.get_is_reverse(c) ? source.flip(c) : c),
                                        source.get_id(c));
                                next_handles.push_back(x);
                                accumulated_length += subgraph.get_length(x);
                            }
                        });
                    }
                }
                curr_handles = std::move(next_handles);
            }
        }

        void extract_id_range(const graph_t &source, const nid_t &id1, const nid_t &id2, graph_t &subgraph,
                              const std::string &progress_message) {
            const bool show_progress = !progress_message.empty();

            std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress;
            if (show_progress) {
                progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(id2 - id1 + 1, progress_message);
            }

            for (nid_t id = id1; id <= id2; ++id) {
                if (!subgraph.has_node(id)) {
                    const handle_t cur_handle = source.get_handle(id);
                    subgraph.create_handle(source.get_sequence(cur_handle), id);
                }

                if (show_progress) {
                    progress->increment(1);
                }
            }

            if (show_progress) {
                progress->finish();
            }
        }

        void embed_lace_paths(graph_t &source, graph_t &subgraph,
                              const std::vector<path_handle_t>& lace_paths) {

            std::vector<std::pair<step_handle_t, step_handle_t>> ranges_to_lace;
            for (auto& path : lace_paths) {
                // check if the nodes are in the output subgraph
                bool in_match = true;
                step_handle_t lace_start = source.path_begin(path);
                // get each range that isn't included
                source.for_each_step_in_path(
                    path,
                    [&](const step_handle_t& step) {
                        if (!subgraph.has_node(source.get_id(source.get_handle_of_step(step)))) {
                            if (in_match) {
                                lace_start = step;
                                in_match = false;
                            }
                        } else {
                            if (!in_match) {
                                ranges_to_lace.push_back(std::make_pair(lace_start, step));
                            }
                            in_match = true;
                        }
                    });
                if (!in_match) {
                    ranges_to_lace.push_back(std::make_pair(lace_start, source.path_back(path)));
                }
            }
            for (auto& range : ranges_to_lace) {
                // get its sequence
                if (range.first != range.second) {
                    std::string seq;
                    for (step_handle_t step = range.first;
                         step != range.second;
                         step = source.get_next_step(step)) {
                        seq += source.get_sequence(source.get_handle_of_step(step));
                    }
                    // add a node with this sequence to both graphs using the same id
                    assert(seq.size());
                    const auto h = source.create_handle(seq);
                    subgraph.create_handle(seq, source.get_id(h));
                    // rewrite the segment in the source graph (nb. inclusive range)
                    source.rewrite_segment(range.first, source.get_previous_step(range.second), { h });
                }
            }
        }

        void windows_in_out(
                const PathHandleGraph& graph,
                const std::vector<path_handle_t>& paths,
                const std::function<bool(handle_t)>& in_bounds,
                const uint64_t& length,
                const std::function<void(const std::vector<path_range_t>&)>& output,
                const uint64_t& num_threads) {

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
            for (path_handle_t path : paths) {
                std::vector<path_range_t> path_ranges;

                uint64_t walked = 0;
                bool first_node = true;

                uint64_t start, end;
                graph.for_each_step_in_path(path, [&](const step_handle_t &step) {
                    const handle_t handle = graph.get_handle_of_step(step);
                    const uint64_t handle_length = graph.get_length(handle);

                    if (in_bounds(handle)) {
                        if (first_node) {
                            first_node = false;
                            start = walked;
                        }

                        end = walked + handle_length;
                    } else if (!first_node) {
                        if (path_ranges.size() > 0 && (start - path_ranges.back().end.offset) < length) {
                            path_ranges.back().end.offset = end;
                        } else {
                            path_ranges.push_back({
                                                          {path, start, false},
                                                          {path, end, false},
                                                          false,
                                                          ""
                                                  });
                        }

                        first_node = true;
                    }

                    walked += handle_length;
                });

                // last path range
                if (!first_node) {
                    if (path_ranges.size() > 0 && (start - path_ranges.back().end.offset) < length) {
                        path_ranges.back().end.offset = end;
                    } else {
                        path_ranges.push_back({
                                                      {path, start, false},
                                                      {path, end, false},
                                                      false,
                                                      ""
                                              });
                    }
                }

                output(path_ranges);
            }
        }

        bool check_and_get_windows_in_out_parameter(
                const std::string& parameter,
                uint64_t &windows_len, uint64_t &windows_min, uint64_t &windows_max, bool &only_tips) {
            const std::regex regex(":");
            const std::vector<std::string> splitted(
                    std::sregex_token_iterator(parameter.begin(), parameter.end(), regex, -1),
                    std::sregex_token_iterator()
            );

            if (splitted.size() != 4) {
                return false;
            }
            if (!utils::is_number(splitted[0])
                || !utils::is_number(splitted[1])
                || !utils::is_number(splitted[2])
                || !utils::is_number(splitted[3])) {
                return false;
            }
            if (stoull(splitted[1]) > stoull(splitted[2])) {
                return false;
            }

            windows_len = stoull(splitted[0]);
            windows_min = stoull(splitted[1]);
            windows_max = stoull(splitted[2]);
            only_tips = (stoull(splitted[3]) == 1);
            return true;
        }
    }
}
