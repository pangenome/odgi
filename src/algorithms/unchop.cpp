/**
 * \file unchop.cpp
 *
 * Defines an algorithm to join adjacent handles.
 */

#include "unchop.hpp"

namespace odgi {
    namespace algorithms {

        /// Concatenates the nodes into a new node with the same external linkage as
        /// the provided component. All handles must be in left to right order and a
        /// consistent orientation. All paths present must run all the way through the
        /// run of nodes from start to end or end to start.
        ///
        /// Returns the handle to the newly created node in the new graph.
        ///
        /// After calling this on a vg::VG, paths will be invalid until
        /// Paths::compact_ranks() is called.
        handle_t concat_nodes(handlegraph::MutablePathDeletableHandleGraph &graph, const std::vector<handle_t> &nodes) {

            // Make sure we have at least 2 nodes
            assert(!nodes.empty() && nodes.front() != nodes.back());

            /*
        #ifdef debug
            std::cerr << "Paths before concat: " << std::endl;
            graph.for_each_path_handle([&](const path_handle_t p) {
               std::cerr << graph.get_path_name(p) << ": ";
               for (auto h : graph.scan_path(p)) {
                   std::cerr << graph.get_id(h) << (graph.get_is_reverse(h) ? '-' : '+') << " ";
               }
               std::cerr << std::endl;
           });
        #endif
            */

            // We also require no edges enter or leave the run of nodes, but we can't check that now.

            // Make the new node
            handle_t new_node;
            {
                std::stringstream ss;
                for (auto &n : nodes) {
                    ss << graph.get_sequence(n);
                }

                new_node = graph.create_handle(ss.str());
            }

#ifdef debug
            std::cerr << "Concatenating ";
            for (auto& n : nodes) {
                std::cerr << graph.get_id(n) << (graph.get_is_reverse(n) ? "-" : "+") << " ";
            }
            std::cerr << "into " << graph.get_id(new_node) << "+" << std::endl;
#endif

            // We should be able to rely on our handle graph to deduplicate edges, but see https://github.com/vgteam/libbdsg/issues/39
            // So we deduplicate ourselves.

            // Find all the neighbors. Make sure to translate edges to the other end of
            // the run, or self loops.
            ska::flat_hash_set<handle_t> left_neighbors;
            graph.follow_edges(nodes.front(), true, [&](const handle_t &left_neighbor) {
                if (left_neighbor == nodes.back()) {
                    // Loop back to the end
                    left_neighbors.insert(new_node);
                } else if (left_neighbor == graph.flip(nodes.front())) {
                    // Loop back to the front
                    left_neighbors.insert(graph.flip(new_node));
                } else {
                    // Normal edge
                    left_neighbors.insert(left_neighbor);
                }
            });

            ska::flat_hash_set<handle_t> right_neighbors;
            graph.follow_edges(nodes.back(), false, [&](const handle_t &right_neighbor) {
                if (right_neighbor == nodes.front()) {
                    // Loop back to the front.
                    // We will have seen it from the other side, so ignore it here.
                } else if (right_neighbor == graph.flip(nodes.back())) {
                    // Loop back to the end
                    right_neighbors.insert(graph.flip(new_node));
                } else {
                    // Normal edge
                    right_neighbors.insert(right_neighbor);
                }
            });

            // Make all the edges, now that we can't interfere with edge listing
            for (auto &n : left_neighbors) {
#ifdef debug
                std::cerr << "Creating edge " << graph.get_id(n) << (graph.get_is_reverse(n) ? "-" : "+") << " -> "
                          <<  graph.get_id(new_node) << (graph.get_is_reverse(new_node) ? "-" : "+") << std::endl;
#endif
                graph.create_edge(n, new_node);
            }
            for (auto &n : right_neighbors) {

#ifdef debug
                std::cerr << "Creating edge " << graph.get_id(new_node) << (graph.get_is_reverse(new_node) ? "-" : "+") << " -> "
                          <<  graph.get_id(n) << (graph.get_is_reverse(n) ? "-" : "+") << std::endl;
#endif

                graph.create_edge(new_node, n);
            }

            {
                // Collect the first and last visits along paths. TODO: this requires
                // walking each path all the way through.
                //
                // This contains the first and last handles in path orientation, and a flag
                // for if the path runs along the reverse strand of our run of nodes.
                std::vector<std::tuple<step_handle_t, step_handle_t, bool>> ranges_to_rewrite;

                graph.for_each_step_on_handle(nodes.front(), [&](const step_handle_t &front_step) {
                    auto path = graph.get_path_handle_of_step(front_step);
#ifdef debug
                    std::cerr << "Consider path " << graph.get_path_name(path) << std::endl;
#endif

                    // If we don't get the same oriented node as where this step is
                    // stepping, we must be stepping on the other orientation.
                    bool runs_reverse = (graph.get_handle_of_step(front_step) != nodes.front());

                    step_handle_t back_step = front_step;

                    while (graph.get_handle_of_step(back_step) !=
                           (runs_reverse ? graph.flip(nodes.back()) : nodes.back())) {
                        // Until we find the step on the path that visits the last node in our run.
                        // Go along the path towards where our last node should be, in our forward orientation.
                        back_step = runs_reverse ? graph.get_previous_step(back_step) : graph.get_next_step(back_step);
                    }

                    // Now we can record the range to rewrite
                    // Make sure to put it into path-forward order
                    if (runs_reverse) {
#ifdef debug
                        std::cerr << "\tGoing to rewrite between " << graph.get_id(graph.get_handle_of_step(front_step)) << " and " << graph.get_id(graph.get_handle_of_step(back_step)) << " backward" << std::endl;
#endif
                        ranges_to_rewrite.emplace_back(back_step, front_step, true);
                    } else {

#ifdef debug
                        std::cerr << "\tGoing to rewrite between " << graph.get_id(graph.get_handle_of_step(front_step)) << " and " << graph.get_id(graph.get_handle_of_step(back_step)) << std::endl;
#endif
                        ranges_to_rewrite.emplace_back(front_step, back_step, false);
                    }
                });

                uint64_t i = 0;
                for (auto &range : ranges_to_rewrite) {
                    // Rewrite each range to visit the new node in the appropriate orientation instead of whatever it did before
                    // Make sure to advance the end of the range because rewrite is end-exclusive (to allow insert).
                    graph.rewrite_segment(std::get<0>(range), std::get<1>(range),
                                          {std::get<2>(range) ? graph.flip(new_node) : new_node});
                }
            }

            // Destroy all the old edges
            // We know they only exist to the left and right neighbors, and along the run
            for (auto &n : left_neighbors) {
                graph.destroy_edge(n, nodes.front());
            }
            for (auto &n : right_neighbors) {
                graph.destroy_edge(nodes.back(), n);
            }
            auto it = nodes.begin();
            auto next_it = it;
            ++next_it;
            while (next_it != nodes.end()) {
                graph.destroy_edge(*it, *next_it);
                it = next_it;
                ++next_it;
            }

            for (auto &n : nodes) {
                // Destroy all the old nodes
#ifdef debug
                std::cerr << "Destroying node " << graph.get_id(n) << std::endl;
#endif
                graph.destroy_handle(n);
            }

            /*
              #ifdef debug
              std::cerr << "Paths after concat: " << std::endl;

              graph.for_each_path_handle([&](const path_handle_t p) {
              std::cerr << graph.get_path_name(p) << ": ";
              for (auto h : graph.scan_path(p)) {
              std::cerr << graph.get_id(h) << (graph.get_is_reverse(h) ? '-' : '+') << " ";
              }
              std::cerr << std::endl;
              });

              #endif
            */

            // Return the new handle we merged to.
            return new_node;
        }

        handle_t combine_handles(handlegraph::MutablePathDeletableHandleGraph &graph,
                                 const std::vector<handle_t> &handles) {
            std::string seq;
            for (auto &handle : handles) {
                seq.append(graph.get_sequence(handle));
            }
            handle_t combined = graph.create_handle(seq);
            // relink the inbound and outbound nodes
            // get the edge context
            std::vector<handle_t> edges_fwd_fwd;
            std::vector<handle_t> edges_fwd_rev;
            std::vector<handle_t> edges_rev_fwd;
            std::vector<handle_t> edges_rev_rev;
            graph.follow_edges(
                    handles.back(), false,
                    [&](const handle_t &h) {
                        edges_fwd_fwd.push_back(h);
                    });
            graph.follow_edges(
                    handles.front(), true,
                    [&](const handle_t &h) {
                        edges_fwd_rev.push_back(h);
                    });
            // destroy the old handles
            for (auto &handle : handles) {
                graph.destroy_handle(handle);
            }
            // connect the ends to the previous context
            // check that we're not trying to make edges that connect back with the nodes in the component
            // there are three cases
            // self looping, front and rear inverting
            for (auto &h : edges_fwd_fwd) {
                if (h == handles.front()) {
                    graph.create_edge(combined, combined);
                } else if (h == graph.flip(handles.back())) {
                    graph.create_edge(combined, graph.flip(combined));
                } else {
                    graph.create_edge(combined, h);
                }
            }
            for (auto &h : edges_fwd_rev) {
                if (h == handles.back()) {
                    graph.create_edge(combined, combined);
                } else if (h == graph.flip(handles.front())) {
                    graph.create_edge(graph.flip(combined), combined);
                } else {
                    graph.create_edge(h, combined);
                }
            }
            return combined;
        }

        bool unchop(handlegraph::MutablePathDeletableHandleGraph &graph) {
            return unchop(graph, 1, false);
        }

        bool unchop(handlegraph::MutablePathDeletableHandleGraph &graph,
                    const uint64_t &nthreads,
                    const bool &show_info) {
#ifdef debug
            std::cerr << "Running unchop" << std::endl;
#endif

            ska::flat_hash_map<nid_t, uint64_t> node_rank;
            uint64_t rank = 0;
            graph.for_each_handle([&](const handle_t &h) {
                node_rank[graph.get_id(h)] = rank++;
            });

            // if possible, don't hold these in memory
            // meanwhile, we exploit the it for parallelizing the paths' validation
            std::vector<std::string> path_names;
            path_names.resize(graph.get_path_count());

            std::vector<std::string> path_sequences;
            path_sequences.resize(graph.get_path_count());

            rank = 0;
            graph.for_each_path_handle([&](const path_handle_t &p) {
                path_names[rank] = graph.get_path_name(p);

                auto &seq = path_sequences[rank];
                graph.for_each_step_in_path(p, [&](const step_handle_t &s) {
                    seq.append(graph.get_sequence(graph.get_handle_of_step(s)));
                });

                ++rank;
            });

            auto components = simple_components(graph, 2, true, nthreads);
            ska::flat_hash_set<nid_t> to_merge;
            for (auto &comp : components) {
                for (auto &handle : comp) {
                    to_merge.insert(graph.get_id(handle));
                }
            }
            std::vector<std::pair<double, handle_t>> ordered_handles;
            graph.for_each_handle(
                    [&](const handle_t &handle) {
                        if (!to_merge.count(graph.get_id(handle))) {
                            ordered_handles.push_back(std::make_pair(
                                    node_rank[graph.get_id(handle)],
                                    handle));
                        }
                    });

            uint64_t num_node_unchopped = 0;
            uint64_t num_new_nodes = 0;
            for (auto &comp : components) {
#ifdef debug
                std::cerr << "Unchop " << comp.size() << " nodes together" << std::endl;
#endif
                if (comp.size() >= 2) {
                    // sort by lowest rank to maintain order
                    double rank_sum = 0;
                    for (auto &handle : comp) {
                        rank_sum += node_rank[graph.get_id(handle)];
                    }
                    double rank_v = rank_sum / comp.size();
                    handle_t n = concat_nodes(graph, comp);
                    ordered_handles.push_back(std::make_pair(rank_v, n));
                    //node_order.push_back(graph.get_id(n));
                    num_node_unchopped += comp.size();
                    num_new_nodes++;
                } else {
                    for (auto &c : comp) {
                        ordered_handles.push_back(std::make_pair(node_rank[graph.get_id(c)], c));
                    }
                }
            }

            // todo try sorting again

            if (show_info) {
                std::cerr << "[odgi::unchop] unchopped " << num_node_unchopped << " nodes into " << num_new_nodes
                          << " new nodes." << std::endl;
            }

            assert(graph.get_node_count() == ordered_handles.size());

            ips4o::parallel::sort(ordered_handles.begin(), ordered_handles.end(), std::less<>(), nthreads);

            std::vector<handle_t> handle_order;
            for (auto &h : ordered_handles) {
                handle_order.push_back(h.second);
            }

            graph.apply_ordering(handle_order, true);

            std::atomic<bool> ok(true);

#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads)
            for (rank = 0; rank < path_names.size(); ++rank) {
                std::string seq;
                graph.for_each_step_in_path(
                        graph.get_path_handle(path_names[rank]),
                        [&](const step_handle_t &s) {
                            seq.append(graph.get_sequence(graph.get_handle_of_step(s)));
                        });

                if (seq != path_sequences[rank]) {
                    std::cerr << "[odgi::algorithms::unchop] failure in unchop" << std::endl;
                    std::cerr << ">expected_" << path_names[rank] << std::endl
                              << path_sequences[rank] << std::endl
                              << ">got_" << path_names[rank] << std::endl << seq << std::endl;

                    ok.store(false);
                }

                std::string().swap(path_names[rank]);
                std::string().swap(path_sequences[rank]);
            }

            return ok.load();
        }
    }
}
