#include "explode.hpp"

namespace odgi {
    namespace algorithms {

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
            add_connecting_edges_to_subgraph(source, subgraph);
        }

        // Create a subpath name
        string Paths_make_subpath_name(const string &path_name, size_t offset, size_t end_offset) {
            string out_name = path_name + "[" + std::to_string(offset);
            if (end_offset > 0) {
                out_name += "-" + std::to_string(end_offset);
            }
            out_name += "]";
            return out_name;
        }

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

    }

}
