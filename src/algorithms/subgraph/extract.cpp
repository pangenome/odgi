#include "extract.hpp"

namespace odgi {
    namespace algorithms {

        std::pair <step_handle_t, int64_t>
        get_graph_pos(const odgi::graph_t &graph, const path_handle_t path, const uint64_t &pos) {
            auto path_end = graph.path_end(path);
            uint64_t walked = 0;
            for (step_handle_t s = graph.path_begin(path);
                 s != path_end; s = graph.get_next_step(s)) {
                handle_t h = graph.get_handle_of_step(s);
                uint64_t node_length = graph.get_length(h);
                if (walked + node_length > pos) {
                    return std::make_pair(s, walked);
                }
                walked += node_length;
            }
#pragma omp critical (cout)
            std::cerr << "[odgi::get_graph_pos] warning: position " << graph.get_path_name(path) << ":" << pos
                      << " outside of path" << std::endl;
            return std::make_pair(path_end, -1);
        }

        // Create a subpath name
        string make_subpath_name(const string &path_name, size_t offset, size_t end_offset) {
            string out_name = path_name + "[" + std::to_string(offset);
            if (end_offset > 0) {
                out_name += "-" + std::to_string(end_offset);
            }
            out_name += "]";
            return out_name;
        }

        void add_subpaths_to_subgraph(const graph_t &source, graph_t &subgraph, bool subpath_naming) {
            auto get_position_of_step =
                    [](const odgi::graph_t &graph, const step_handle_t &step_to_find) {
                        path_handle_t path = graph.get_path_handle_of_step(step_to_find);
                        auto path_end = graph.path_end(path);
                        uint64_t walked = 0;
                        for (step_handle_t s = graph.path_begin(path); s != path_end; s = graph.get_next_step(s)) {
                            if (s == step_to_find) {
                                return walked;
                            }

                            handle_t h = graph.get_handle_of_step(s);
                            uint64_t node_length = graph.get_length(h);
                            walked += node_length;
                        }
#pragma omp critical (cout)
                        std::cerr << "[odgi::get_position_of_step] warning: step "
                                  << graph.get_id((graph.get_handle_of_step(step_to_find))) << " in "
                                  << graph.get_path_name(path) << " not found" << std::endl;

                        return walked;
                    };

            std::unordered_map <std::string, std::map<uint64_t, handle_t>> subpaths;
            subgraph.for_each_handle([&](const handle_t &h) {
                handlegraph::nid_t id = subgraph.get_id(h);
                if (source.has_node(id)) {
                    handle_t handle = source.get_handle(id);
                    source.for_each_step_on_handle(handle, [&](const step_handle_t &step) {
                        path_handle_t path = source.get_path_handle_of_step(step);
                        std::string path_name = source.get_path_name(path);
                        uint64_t pos = get_position_of_step(source, step);
                        subpaths[path_name][pos] = source.get_is_reverse(source.get_handle_of_step(step))
                                                   ? subgraph.flip(h) : h;

                    });
                    /*//return iteratee(step, get_is_reverse(get_handle_of_step(step)) != get_is_reverse(handle), get_position_of_step(step));
                    source.for_each_step_position_on_handle(handle, [&](const step_handle_t& step, const bool& is_rev, const uint64_t& pos) {
                        path_handle_t path = source.get_path_handle_of_step(step);
                        std::string path_name = source.get_path_name(path);
                        subpaths[path_name][pos] = is_rev ? subgraph.flip(h) : h;
                        return true;
                    });*/
                }
            });

            function < path_handle_t(
            const string&, bool, size_t)> new_subpath = [&subgraph](const string &path_name, bool is_circular,
                                                                    size_t subpath_offset) {
                string subpath_name = make_subpath_name(path_name, subpath_offset,
                                                        0);
                if (subgraph.has_path(subpath_name)) {
                    subgraph.destroy_path(subgraph.get_path_handle(subpath_name));
                }
                return subgraph.create_path_handle(subpath_name, is_circular);
            };

            for (auto &subpath : subpaths) {
                const std::string &path_name = subpath.first;
                path_handle_t source_path_handle = source.get_path_handle(path_name);
                // destroy the path if it exists
                if (subgraph.has_path(path_name)) {
                    subgraph.destroy_path(subgraph.get_path_handle(path_name));
                }
                // create a new path.  give it a subpath name if the flag's on and its smaller than original
                path_handle_t path;
                if (!subpath_naming || subpath.second.size() == source.get_step_count(source_path_handle) ||
                    subpath.second.empty()) {
                    path = subgraph.create_path_handle(path_name, source.get_is_circular(source_path_handle));
                } else {
                    path = new_subpath(path_name, source.get_is_circular(source_path_handle),
                                       subpath.second.begin()->first);
                }
                for (auto p = subpath.second.begin(); p != subpath.second.end(); ++p) {
                    const handle_t &handle = p->second;
                    if (p != subpath.second.begin() && subpath_naming) {
                        auto prev = p;
                        --prev;
                        const handle_t &prev_handle = prev->second;
                        // distance from map
                        size_t delta = p->first - prev->first;
                        // what the distance should be if they're contiguous depends on relative orienations
                        size_t cont_delta = subgraph.get_length(prev_handle);
                        if (delta != cont_delta) {
                            // we have a discontinuity!  we'll make a new path can continue from there
                            assert(subgraph.get_step_count(path) > 0);
                            path = new_subpath(path_name, subgraph.get_is_circular(path), p->first);
                        }
                    }
                    //fill in the path information
                    subgraph.append_step(path, handle);
                }
            }
        }

        void extract_path_range(const graph_t &source, path_handle_t path_handle, int64_t start, int64_t end,
                                graph_t &subgraph) {

            // ToDo: could be improved getting and extracting simultaneously (not first get_graph_pos and then extract_path_range)?
            std::pair <step_handle_t, uint64_t> start_step_and_start_pos = get_graph_pos(source, path_handle, start);

            step_handle_t start_step = start_step_and_start_pos.first;  //source.get_step_at_position(path_handle, start);
            size_t start_position = start_step_and_start_pos.second;    //source.get_position_of_step(start_step);
            size_t size_needed = end < 0 ? numeric_limits<size_t>::max() : end - start + 1 + start - start_position;
            size_t running_length = 0;

            auto path_end = source.path_end(path_handle);
            for (step_handle_t cur_step = start_step;
                 cur_step != path_end && running_length < size_needed; cur_step = source.get_next_step(cur_step)) {

                handle_t cur_handle = source.get_handle_of_step(cur_step);
                subgraph.create_handle(source.get_sequence(cur_handle), source.get_id(cur_handle));
                if (cur_step != start_step) {
                    handle_t prev_handle = source.get_handle_of_step(source.get_previous_step(cur_step));
                    subgraph.create_edge(
                            subgraph.get_handle(source.get_id(prev_handle), source.get_is_reverse(prev_handle)),
                            subgraph.get_handle(source.get_id(cur_handle), source.get_is_reverse(cur_handle)));
                }
                running_length += source.get_length(cur_handle);

            }
        }


        /// We can accumulate a subgraph without accumulating all the edges between its nodes
        /// this helper ensures that we get the full set
        void add_connecting_edges_to_subgraph(const graph_t &source, graph_t &subgraph) {
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
            });
        }


        void expand_subgraph_by_steps(const graph_t &source, graph_t &subgraph, const uint64_t &steps,
                                      bool forward_only) {
            std::vector <handle_t> curr_handles;
            subgraph.for_each_handle([&](const handle_t &h) {
                curr_handles.push_back(h);
            });
            for (uint64_t i = 0; i < steps && !curr_handles.empty(); ++i) {
                std::vector <handle_t> next_handles;
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
            add_connecting_edges_to_subgraph(source, subgraph);
        }


        void add_full_paths_to_component(const graph_t &source, graph_t &component) {

            // We want to track the path names in each component
            set <path_handle_t> paths;

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

    }
}
