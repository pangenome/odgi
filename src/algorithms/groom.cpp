/**
 * \file groom.cpp
 *
 * Defines an algorithm to remove spurious inverting links from the graph
 * by exploring the graph from the orientation supported by the most paths.
 */

#include "groom.hpp"

namespace odgi {
    namespace algorithms {

    std::vector<handle_t>groom(const handlegraph::MutablePathDeletableHandleGraph &graph,
							   bool progress_reporting, const std::vector<handlegraph::path_handle_t> target_paths, bool use_bfs) {
			bool target_grooming = (target_paths.size() > 0);

            // This (s) is our set of oriented nodes.
            //dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,256,16> > s;
            uint64_t min_handle_rank = std::numeric_limits<uint64_t>::max();
            uint64_t max_handle_rank = 0;
            graph.for_each_handle(
                    [&](const handle_t &found) {
                        uint64_t handle_rank = number_bool_packing::unpack_number(found);
                        min_handle_rank = std::min(min_handle_rank, handle_rank);
                        max_handle_rank = std::max(max_handle_rank, handle_rank);
                    });

            // Start with the heads of the graph.
            // We could also just use the first node of the graph.
            std::vector<handle_t> seeds;
			if (!target_grooming) {
				bool use_heads = true;
				bool use_tails = false;
				if (use_heads) {
					seeds = head_nodes(&graph);
				} else if (use_tails) {
					seeds = tail_nodes(&graph);
				} else {
					handle_t min_handle = number_bool_packing::pack(min_handle_rank, false);
					seeds = {min_handle};
				}
			}

			// do we have target paths which we want to force to have a forward orientation?
			std::vector<bool> is_ref;
			std::vector<bool> needs_flipping;
			if (target_grooming) {
				std::fill_n(std::back_inserter(is_ref), graph.get_node_count(), false);
				std::fill_n(std::back_inserter(needs_flipping), graph.get_node_count(), false);
				std::unique_ptr<progress_meter::ProgressMeter> target_paths_progress;
				if (progress_reporting) {
					std::string banner = "[odgi::groom] preparing target path vectors:";
					target_paths_progress = std::make_unique<progress_meter::ProgressMeter>(target_paths.size(), banner);
				}
				for (handlegraph::path_handle_t target_path : target_paths) {
					graph.for_each_step_in_path(
							target_path,
							[&](const step_handle_t& step) {
								handle_t handle = graph.get_handle_of_step(step);
								uint64_t i = number_bool_packing::unpack_number(handle);
								if (!is_ref[i]) {
									is_ref[i] = true;
									seeds.push_back(handle);
									// do we need flipping?
									if (graph.get_is_reverse(handle)) {
										needs_flipping[i] = true;
									}
								}
							});
					if (progress_reporting) {
						target_paths_progress->increment(1);
					}
				}
				if (progress_reporting) {
					target_paths_progress->finish();
				}
			}

            // We need to keep track of the nodes we haven't visited to seed subsequent runs of the BFS
            dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector, 256, 16> > unvisited, flipped;
            for (uint64_t i = 0; i <= max_handle_rank; ++i) {
                unvisited.push_back(1);
                flipped.push_back(0);
            }

            uint64_t prev_max_root = 0;
            uint64_t prev_max_length = 0;

            std::unique_ptr<progress_meter::ProgressMeter> bfs_progress;
            if (progress_reporting) {
                std::string banner = "[odgi::groom] grooming:";
                bfs_progress = std::make_unique<progress_meter::ProgressMeter>(graph.get_node_count(), banner);
            }

            uint64_t edge_count = 0;

            while (unvisited.rank1(unvisited.size()) != 0) {
                if (use_bfs) {
                    bfs(graph,
                        [&graph, &unvisited, &flipped, &progress_reporting, &bfs_progress, &needs_flipping, &is_ref, &target_grooming]
                        (const handle_t &h, const uint64_t &r, const uint64_t &l, const uint64_t &d) {
                            if (progress_reporting) {
                                bfs_progress->increment(1);
                            }
                            uint64_t i = number_bool_packing::unpack_number(h);
                            unvisited.set(i, 0);
							if (target_grooming && is_ref[i]) {
								if (needs_flipping[i]) {
									flipped.set(i, true);
								} else {
									flipped.set(i, false);
								}
							} else {
								flipped.set(i, graph.get_is_reverse(h));
							}
                        },
                        [&unvisited](const handle_t &h) {
                            uint64_t i = number_bool_packing::unpack_number(h);
                            return unvisited.at(i) == 0;
                        },
                        [&edge_count](const handle_t &l, const handle_t &h) {
                            ++edge_count;
                            return false;
                        },
                        [](void) { return false; },
                        seeds,
                        {},
                        false); // don't use bidirectional search
                    // get another seed
                    if (unvisited.rank1(unvisited.size()) != 0) {
                        uint64_t i = unvisited.select1(0);
                        handle_t h = number_bool_packing::pack(i, false);
                        seeds = {h};
                    }
                } else {
                    dfs(graph,
                        [&graph, &unvisited, &flipped, &progress_reporting, &bfs_progress, &needs_flipping, &is_ref, &target_grooming]
                        (const handle_t &h) {
                            if (progress_reporting) {
                                bfs_progress->increment(1);
                            }
                            uint64_t i = number_bool_packing::unpack_number(h);
                            unvisited.set(i, 0);
							if (target_grooming && is_ref[i]) {
								if (needs_flipping[i]) {
									flipped.set(i, true);
								} else {
									flipped.set(i, false);
								}
							} else {
								flipped.set(i, graph.get_is_reverse(h));
							}
                        },
                        [&unvisited](const handle_t &h) {
                            uint64_t i = number_bool_packing::unpack_number(h);
                            return unvisited.at(i) == 0;
                        },
                        [](const handle_t& h) { return false; },
                        [](void) { return false; },
                        seeds);
                    //{},
                    //false); // don't use bidirectional search
                    // get another seed
                    if (unvisited.rank1(unvisited.size()) != 0) {
                        uint64_t i = unvisited.select1(0);
                        handle_t h = number_bool_packing::pack(i, false);
                        seeds = {h};
                    }
                }
            }

            if (progress_reporting) {
                bfs_progress->finish();
            }

            std::unique_ptr<progress_meter::ProgressMeter> handle_progress;
            if (progress_reporting) {
                std::string banner = "[odgi::groom] organizing handles:";
                handle_progress = std::make_unique<progress_meter::ProgressMeter>(graph.get_node_count(), banner);
            }

            std::vector<handle_t> order;
            uint64_t num_flipped_handles = 0;

            graph.for_each_handle(
                    [&graph, &order, &flipped, &num_flipped_handles, &progress_reporting, &handle_progress](
                            const handle_t &h) {
                        bool to_flip = flipped[number_bool_packing::unpack_number(h)];
                        if (!to_flip) {
                            order.push_back(h);
                        } else {
                            order.push_back(graph.flip(h));
                            ++num_flipped_handles;
                        }

                        if (progress_reporting) {
                            handle_progress->increment(1);
                        }
                    });

            if (progress_reporting) {
                handle_progress->finish();

                std::cerr << "[odgi::groom] flipped " << num_flipped_handles << " handles" << std::endl;
            }

            return order;

        }

    }
}
