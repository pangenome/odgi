#include "extract.hpp"

namespace odgi {
    namespace algorithms {

        void add_full_paths_to_component(const graph_t &source, graph_t &component) {

            // We want to track the path names in each component
            set<path_handle_t> paths;

            // Record paths
            component.for_each_handle([&](const handle_t &h) {
                handlegraph::nid_t id = component.get_id(h);

                if (source.has_node(id)) {
                    handle_t source_handle = source.get_handle(id);

                    source.for_each_step_on_handle(source_handle, [&](const step_handle_t &source_step) {
                        paths.insert(source.get_path_handle_of_step(source_step));
                    });
                }
            });

            // Copy the paths over
            for (path_handle_t path_handle : paths) {
                path_handle_t new_path_handle = component.create_path_handle(source.get_path_name(path_handle),
                                                                             source.get_is_circular(path_handle));
                for (handle_t handle : source.scan_path(path_handle)) {
                    component.append_step(new_path_handle, component.get_handle(source.get_id(handle),
                                                                                source.get_is_reverse(handle)));
                }
            }
        }

        // Create a subpath name
        string make_subpath_name(const string &path_name, size_t offset, size_t end_offset) {
            string out_name = path_name + ":" + std::to_string(offset);
            if (end_offset > offset) {
                out_name += "-" + std::to_string(end_offset);
            }
            return out_name;
        }

        void add_subpaths_to_subgraph(const graph_t &source, const std::vector<path_handle_t> source_paths, graph_t &subgraph, const uint64_t num_threads, const std::string& progress_message) {
            bool show_progress = !progress_message.empty();

            std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress;
            if (show_progress) {
                progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                        source_paths.size(), progress_message);
            }

            auto create_and_fill_subpath = [](graph_t &subgraph, const string &path_name,
                                              const bool is_circular, const size_t start, const size_t end,
                                              std::vector<handle_t>& handles_to_embed) {
                string subpath_name = make_subpath_name(path_name, start, end);
                if (subgraph.has_path(subpath_name)) {
                    subgraph.destroy_path(subgraph.get_path_handle(subpath_name));
                }
                path_handle_t path = subgraph.create_path_handle(subpath_name, is_circular);

                for (auto h_to_embed : handles_to_embed) {
                    subgraph.append_step(path, h_to_embed);
                }

                std::vector<handle_t>().swap(handles_to_embed);
            };

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
            for (auto source_path_handle : source_paths) {
                std::string path_name = source.get_path_name(source_path_handle);

                uint64_t walked = 0;

                uint64_t start, end;
                std::vector<handle_t> handles_to_embed;
                source.for_each_step_in_path(source_path_handle, [&](const step_handle_t &step) {
                    handle_t source_handle = source.get_handle_of_step(step);
                    uint64_t source_length = source.get_length(source_handle);
                    handlegraph::nid_t source_id = source.get_id(source_handle);

                    if (subgraph.has_node(source_id)) {
                        if (handles_to_embed.empty()) {
                            start = walked;
                        }

                        end = walked + source_length;
                        handles_to_embed.push_back(
                                subgraph.get_handle(source_id, source.get_is_reverse(source_handle)));
                    } else if (!handles_to_embed.empty()) {
                        create_and_fill_subpath(subgraph, path_name, source.get_is_circular(source_path_handle), start,
                                                end, handles_to_embed);
                    }

                    walked += source_length;
                });

                // last subpath
                if (!handles_to_embed.empty()) {
                    create_and_fill_subpath(subgraph, path_name, source.get_is_circular(source_path_handle), start, end,
                                            handles_to_embed);
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
            //bool first_step = true;
            uint64_t walked = 0;
            auto path_end = source.path_end(path_handle);
            for (step_handle_t cur_step = source.path_begin(path_handle);
                 cur_step != path_end && walked <= end; cur_step = source.get_next_step(cur_step)) {
                handle_t cur_handle = source.get_handle_of_step(cur_step);
                walked += source.get_length(cur_handle);
                if (walked > start) {
                    nid_t id = source.get_id(cur_handle);
                    if (!subgraph.has_node(id)) {
                        subgraph.create_handle(source.get_sequence(source.get_handle(id)), id);
                    }

                    // This code ensures that there are edges connecting the nodes the path crosses.
                    // Commenting it, it is assumed that the topology of the graph is consistent with the path
                    /*if (!first_step) {
                        handle_t prev_handle = source.get_handle_of_step(source.get_previous_step(cur_step));
                        subgraph.create_edge(
                                subgraph.get_handle(source.get_id(prev_handle), source.get_is_reverse(prev_handle)),
                                subgraph.get_handle(source.get_id(cur_handle), source.get_is_reverse(cur_handle)));
                    } else {
                        first_step = false;
                    }*/
                }
            }
        }

        /// We can accumulate a subgraph without accumulating all the edges between its nodes
        /// this helper ensures that we get the full set
        void add_connecting_edges_to_subgraph(const graph_t &source, graph_t &subgraph, const std::string& progress_message) {
            bool show_progress = !progress_message.empty();

            std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress;
            if (show_progress) {
                progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                        subgraph.get_node_count(), progress_message);
            }

            subgraph.for_each_handle([&](const handle_t &handle) {
                nid_t id = subgraph.get_id(handle);
                handle_t source_handle = source.get_handle(id, subgraph.get_is_reverse(handle));
                source.follow_edges(source_handle, false, [&](const handle_t &next) {
                    nid_t next_id = source.get_id(next);
                    if (subgraph.has_node(next_id)) {
                        handle_t subgraph_next = subgraph.get_handle(next_id, source.get_is_reverse(next));
                        if (!subgraph.has_edge(handle, subgraph_next)) {
                            subgraph.create_edge(handle, subgraph_next);
                        }
                    }
                });
                source.follow_edges(source_handle, true, [&](const handle_t &prev) {
                    nid_t prev_id = source.get_id(prev);
                    if (subgraph.has_node(prev_id)) {
                        handle_t subgraph_prev = subgraph.get_handle(prev_id, source.get_is_reverse(prev));
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

        void expand_subgraph_by_steps(const graph_t &source, graph_t &subgraph, const uint64_t &steps,
                                      bool forward_only, const std::string& progress_message) {
            bool show_progress = !progress_message.empty();

            std::vector<handle_t> curr_handles;
            subgraph.for_each_handle([&](const handle_t &h) {
                curr_handles.push_back(h);
            });
            for (uint64_t i = 0; i < steps && !curr_handles.empty(); ++i) {
                std::vector<handle_t> next_handles;
                for (auto &h : curr_handles) {
                    handle_t old_h = source.get_handle(subgraph.get_id(h));
                    source.follow_edges(old_h, false, [&](const handle_t &c) {
                        handle_t x;
                        if (!subgraph.has_node(source.get_id(c))) {
                            x = subgraph.create_handle(
                                    source.get_sequence(source.get_is_reverse(c) ? source.flip(c) : c),
                                    source.get_id(c));
                            next_handles.push_back(x);
                        } else {
                            x = subgraph.get_handle(source.get_id(c));
                        }
                        if (source.get_is_reverse(c)) {
                            x = subgraph.flip(x);
                        }
                        if (!subgraph.has_edge(h, x)) {
                            subgraph.create_edge(h, x);
                        }
                    });
                    if (!forward_only) {
                        source.follow_edges(old_h, true, [&](const handle_t &c) {
                            handle_t x;
                            if (!subgraph.has_node(source.get_id(c))) {
                                x = subgraph.create_handle(
                                        source.get_sequence(source.get_is_reverse(c) ? source.flip(c) : c),
                                        source.get_id(c));
                                next_handles.push_back(x);
                            } else {
                                x = subgraph.get_handle(source.get_id(c));
                            }
                            if (source.get_is_reverse(c)) {
                                x = subgraph.flip(x);
                            }
                            if (!subgraph.has_edge(x, h)) {
                                subgraph.create_edge(x, h);
                            }
                        });
                    }
                }
                curr_handles = std::move(next_handles);
            }
            add_connecting_edges_to_subgraph(source, subgraph, progress_message);
        }

        void expand_subgraph_by_length(const graph_t &source, graph_t &subgraph, const uint64_t &length,
                                       bool forward_only, const std::string& progress_message) {
            bool show_progress = !progress_message.empty();

            uint64_t accumulated_length = 0;
            std::vector<handle_t> curr_handles;
            subgraph.for_each_handle([&](const handle_t &h) {
                curr_handles.push_back(h);
            });
            while (accumulated_length < length && !curr_handles.empty()) {
                std::vector<handle_t> next_handles;
                for (auto &h : curr_handles) {
                    handle_t old_h = source.get_handle(subgraph.get_id(h));
                    source.follow_edges(old_h, false, [&](const handle_t &c) {
                        handle_t x;
                        if (!subgraph.has_node(source.get_id(c))) {
                            x = subgraph.create_handle(
                                    source.get_sequence(source.get_is_reverse(c) ? source.flip(c) : c),
                                    source.get_id(c));
                            next_handles.push_back(x);
                            accumulated_length += subgraph.get_length(x);
                        } else {
                            x = subgraph.get_handle(source.get_id(c));
                        }
                        if (source.get_is_reverse(c)) {
                            x = subgraph.flip(x);
                        }
                        if (!subgraph.has_edge(h, x)) {
                            subgraph.create_edge(h, x);
                        }
                    });
                    if (!forward_only) {
                        source.follow_edges(old_h, true, [&](const handle_t &c) {
                            handle_t x;
                            if (!subgraph.has_node(source.get_id(c))) {
                                x = subgraph.create_handle(
                                        source.get_sequence(source.get_is_reverse(c) ? source.flip(c) : c),
                                        source.get_id(c));
                                next_handles.push_back(x);
                                accumulated_length += subgraph.get_length(x);
                            } else {
                                x = subgraph.get_handle(source.get_id(c));
                            }
                            if (source.get_is_reverse(c)) {
                                x = subgraph.flip(x);
                            }
                            if (!subgraph.has_edge(x, h)) {
                                subgraph.create_edge(x, h);
                            }
                        });
                    }
                }
                curr_handles = std::move(next_handles);
            }
            add_connecting_edges_to_subgraph(source, subgraph, progress_message);
        }

        void extract_id_range(const graph_t &source, const nid_t &id1, const nid_t &id2, graph_t &subgraph, const std::string& progress_message) {
            bool show_progress = !progress_message.empty();

            std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress;
            if (show_progress) {
                progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(id2 - id1 + 1, progress_message);
            }

            for (nid_t i = id1; i <= id2; ++i) {
                if (!subgraph.has_node(i)) {
                    subgraph.create_handle(source.get_sequence(source.get_handle(i)), i);
                }

                if (show_progress) {
                    progress->increment(1);
                }
            }

            if (show_progress) {
                progress->finish();
            }
        }
    }
}
