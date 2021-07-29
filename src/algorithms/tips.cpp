#include "tips.hpp"

namespace odgi {
	namespace algorithms {

		using namespace handlegraph;

		void walk_tips(const graph_t& graph,
					   const std::vector<path_handle_t>& paths,
					   const path_handle_t& target_path_t,
					   const std::vector<bool>& target_handles,
					   algorithms::step_index_t& step_index,
					   const uint64_t& num_threads,
					   const std::function<step_handle_t(const path_handle_t&)>& get_path_end,
					   const std::function<step_handle_t(const step_handle_t&)>& get_step,
					   const std::function<bool(const step_handle_t&)>& has_step,
					   algorithms::tips_bed_writer& bed_writer_thread,
					   const bool& progress,
					   const bool& walk_from_front,
					   ska::flat_hash_set<std::string>& not_visited_set,
					   const uint64_t& n_best_mappings,
					   const uint64_t& walking_dist) {

			const std::string target_path = graph.get_path_name(target_path_t);

			std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_meter;
			if (progress) {
				std::string end;
				if (walk_from_front) {
					end = "front";
				} else {
					end = "back";
				}
				std::string progress_message = "[odgi::tips::walk_tips] BED Progress for query path '"
											   + graph.get_path_name(target_path_t) + "' walking from the " + end + " of all paths:";
				progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
						paths.size(), progress_message);
			}

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
			for (auto path : paths) {
				// prevent self tips
				if (path == target_path_t) {
					continue;
				}
				std::string query_path_name = graph.get_path_name(path);
				bool tip_reached_target = false;
				// collecting tips
				step_handle_t cur_step = get_path_end(path);
				handle_t cur_h = graph.get_handle_of_step(cur_step);
				while (!tip_reached_target) {
					// did we already hit the given reference path?
					if (target_handles[number_bool_packing::unpack_number(cur_h)]) {
						// TODO remove this
						std::vector<uint64_t> target_path_start_positions;
						std::vector<uint64_t> target_path_end_positions;
						std::vector<step_handle_t> target_step_handles;
						graph.for_each_step_on_handle(
								cur_h,
								[&](const step_handle_t& s) {
									/// we can do these expensive iterations here, because we only have to do it once for each walk
									if (graph.get_path_handle_of_step(s) == target_path_t) {
										// TODO collect only the steps for the given target
										// we collect all positions
										// later we will get the smallest, the highest and the median from these
										// TODO remove this
										uint64_t pos = step_index.get_position(s);
										target_path_start_positions.push_back(pos);
										target_path_end_positions.push_back(pos + graph.get_length(cur_h) - 1); // 0-based
										target_step_handles.push_back(s);
									}
								});
						/// for the query:
						// walk the given walking_dist first to the left (we might not be able to do the full walk)
						// then walk the given walking_dist to the right
						// in order to collect all the visited nodes and how often they were visited in a ska::flat_hash_map<uint64_t, uint64_t>
						// the key is the node identifier, the value is the number of times this node was visited
						ska::flat_hash_map<nid_t, uint64_t> query_set = collect_nodes_in_walking_dist(graph, walking_dist, cur_step);
//						for (auto pair : query_set) {
//#pragma omp critical (cout)
//							std::cerr << "first: " << pair.first << " second: " << pair.second << std::endl;
//						}
						for (step_handle_t target_step : target_step_handles) {
							ska::flat_hash_map<nid_t, uint64_t> target_set = collect_nodes_in_walking_dist(graph, walking_dist, target_step);
							ska::flat_hash_map<nid_t, uint64_t> union_set;
							for (auto query_item : query_set) {
								union_set[query_item.first] = query_item.second;
							}
							// TODO check sizes of union_set and query_set, hopefully we did not alter the query_set
							std::cerr << "UNION SIZE BEFORE: " << union_set.size() << std::endl;
							add_target_set_to_union_set(union_set, target_set);
							std::cerr << "UNION SIZE AFTER: " << union_set.size() << std::endl;
							/*
							for (auto u : union_set) {
								std::cerr << "node_id: " << u.first << " node_count: " << u.second << std::endl;
							}
							 */
							ska::flat_hash_map<nid_t, uint64_t> intersection_set = intersect_target_query_sets(union_set, target_set, query_set);
							std::cerr << "INTERSECT SIZE: " << intersection_set.size() << std::endl;
							/*
							for (auto i : intersection_set) {
								std::cerr << "node_id: " << i.first << " node_count: " << i.second << std::endl;
							}
							 */
							double jaccard = jaccard_idx_from_intersect_union_sets(intersection_set, union_set, graph);
#pragma omp critical (cout)
							std::cerr << "Jaccard index of query " << query_path_name << " and target " << target_path << ": " << jaccard << ". Direction: " << walk_from_front << std::endl;
						}

//#pragma omp critical (cout)
//					std::cerr << "target_path_start_positions.size(): " << target_path_start_positions.size() << std::endl;
						std::sort(target_path_start_positions.begin(),
								  target_path_start_positions.end(),
								  [&](const uint64_t & pos_a,
									  const uint64_t & pos_b) {
									  return pos_a < pos_b;
								  });
						std::sort(target_path_end_positions.begin(),
								  target_path_end_positions.end(),
								  [&](const uint64_t & pos_a,
									  const uint64_t & pos_b) {
									  return pos_a < pos_b;
								  });
						uint64_t target_min_pos = target_path_start_positions[0]; // 0-based starting position in BED
						uint64_t target_max_pos = target_path_end_positions[target_path_end_positions.size() - 1] + 1; // 1-based ending position in BED
//#pragma omp critical (cout)
						//std::cout << target_path << "\t" << target_min_pos << "\t" << target_max_pos << "\t"
						//		  << target_pos_median << "\t" << graph.get_path_name(path) << "\t" << step_index.get_position(cur_step) << std::endl;
						/// add BED record to queue of BED_writer_thread
						bed_writer_thread.append(target_path, target_min_pos, target_max_pos,
							   query_path_name, step_index.get_position(cur_step), 1.0, walk_from_front);
						tip_reached_target = true;
					}
					if (has_step(cur_step)) {
						cur_step = get_step(cur_step);
						cur_h = graph.get_handle_of_step(cur_step);
					} else {
						// did we iterate over all steps and we did not hit the query path?
						tip_reached_target = true;
#pragma omp critical (not_visited_set)
						not_visited_set.insert(query_path_name);
						// do we still want this?
//						if (progress) {
//#pragma omp critical (cout)
//							std::cerr << "[odgi::tips::walk_tips] warning: For target path '" << query_path_name << "' there was no hit!" << std::endl;
//						}
					}
				}
				if (progress) {
					progress_meter->increment(1);
				}
			}
			if (progress) {
				progress_meter->finish();
			}
		}

		ska::flat_hash_map<nid_t , uint64_t> collect_nodes_in_walking_dist(const graph_t& graph,
																			  const uint64_t& walking_dist,
																			  const step_handle_t& start_step) {
			ska::flat_hash_map<nid_t, uint64_t> node_count_set;
			/// first walk to previous steps up to the walking_dist
			uint64_t dist_walked = 0;
			// where does our step come from
			handle_t cur_h = graph.get_handle_of_step(start_step);
			const bool is_rev = number_bool_packing::unpack_bit(cur_h);
			step_handle_t cur_step = start_step;
			nid_t cur_id = graph.get_id(cur_h);
			if (is_rev) {
				dist_walked += graph.get_length(cur_h);
				node_count_set[cur_id] = 1;
			}
			while (graph.has_previous_step(cur_step) && (dist_walked <= walking_dist)) {
				step_handle_t prev_step = graph.get_previous_step(cur_step);
				handle_t prev_h = graph.get_handle_of_step(prev_step);
				nid_t prev_id = graph.get_id(prev_h);
				//
				if (node_count_set.count(prev_id) == 0) {
					node_count_set[prev_id] = 1;
				} else {
					node_count_set[prev_id] = node_count_set[prev_id] + 1;
				}
				dist_walked += graph.get_length(prev_h);
				cur_step = prev_step;
			}
			dist_walked = 0;
			// where does our step come from
			if (!is_rev) {
				dist_walked += graph.get_length(cur_h);
				if (node_count_set.count(cur_id) == 0) {
					node_count_set[cur_id] = 1;
				} else {
					node_count_set[cur_id] = node_count_set[cur_id] + 1;
				}
			}
			/// walking the next steps up to the walking_dist
			cur_step = start_step;
			while (graph.has_next_step(cur_step) && dist_walked <= walking_dist) {
				step_handle_t next_step = graph.get_next_step(cur_step);
				handle_t next_h = graph.get_handle_of_step(next_step);
				nid_t next_id = graph.get_id(next_h);
				if (node_count_set.count(next_id) == 0) {
					node_count_set[next_id] = 1;
				} else {
					node_count_set[next_id] = node_count_set[next_id] + 1;
				}
				dist_walked += graph.get_length(next_h);
				cur_step = next_step;
			}
			return node_count_set;
		}

		void add_target_set_to_union_set(ska::flat_hash_map<nid_t , uint64_t>& union_set,
								   const ska::flat_hash_map<nid_t , uint64_t>& target_set) {
			for (auto target_item : target_set) {
				nid_t node_id = target_item.first;
				uint64_t node_count = target_item.second;
				if (union_set.count(node_id) != 0) {
					uint64_t union_node_count = union_set[node_id];
					union_set[node_id] = std::max(node_count, union_node_count);
				} else {
					union_set[node_id] = node_count;
				}
			}
		}

		ska::flat_hash_map<nid_t, uint64_t> intersect_target_query_sets(ska::flat_hash_map<nid_t , uint64_t>& union_set,
										 ska::flat_hash_map<nid_t , uint64_t>& target_set,
										 ska::flat_hash_map<nid_t , uint64_t>& query_set) {
			ska::flat_hash_map<nid_t, uint64_t> intersect_set;
			for (auto union_item : union_set) {
				nid_t union_node_id = union_item.first;
				// node is present in both sets, we take the minimum count
				if ((target_set.count(union_node_id) != 0) && (query_set.count(union_node_id) != 0)) {
					uint64_t target_node_count = target_set[union_node_id];
					uint64_t query_node_count = query_set[union_node_id];
					intersect_set[union_node_id] = std::min(target_node_count, query_node_count);
				}
			}
			return intersect_set;
		}

		double jaccard_idx_from_intersect_union_sets(ska::flat_hash_map<nid_t , uint64_t>& intersection_set,
													 ska::flat_hash_map<nid_t , uint64_t>& union_set,
													 const graph_t& graph) {
			double jaccard_idx = 0.0;
			uint64_t intersect_seq_len = 0;
			uint64_t union_seq_len = 0;
			for (auto intersect_item : intersection_set) {
				nid_t intersect_node_id = intersect_item.first;
				uint64_t intersect_node_count = intersect_item.second;
				handle_t h = graph.get_handle(intersect_node_id);
				intersect_seq_len += graph.get_length(h);
			}
			for (auto union_item : union_set) {
				nid_t union_node_id = union_item.first;
				uint64_t union_node_count = union_item.second;
				handle_t h = graph.get_handle(union_node_id);
				union_seq_len += graph.get_length(h);
			}
			std::cerr << "intersect_seq_len: " << intersect_seq_len << std::endl;
			std::cerr << "union_seq_len: " << union_seq_len << std::endl;
			jaccard_idx = (double) intersect_seq_len / (double) union_seq_len;
			return jaccard_idx;
		}

	}
}
