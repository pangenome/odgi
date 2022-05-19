#include "path_jaccard.hpp"

namespace odgi {
	namespace algorithms {

		using namespace handlegraph;

		/// calculate all jaccard indices from a given target_step_handles and a current query step
		std::vector<step_jaccard_t> jaccard_indices_from_step_handles(const graph_t& graph,
																	  const uint64_t& walking_dist,
																	  const step_handle_t& cur_step,
																	  std::vector<step_handle_t>& target_step_handles) {
			/// collect the visited nodes in a vector, so we can do the additional interations faster
			// this is very RAM intensive, maybe it should be made optional
			// gives a speed up of ~1.5x
			// ska::flat_hash_map<step_handle_t, std::pair<std::vector<nid_t>, std::vector<nid_t>>> steps_nodes_prev_next_map;

			/// collect the possible minimal and maximal walking distance in nucleotides from the query step and all possible target steps
			/// we don't care about orientation here
			std::pair<uint64_t , uint64_t> min_max_walk_dist = find_min_max_walk_dist_from_query_targets(graph,
																										 walking_dist,
																										 cur_step,
																										 target_step_handles);
#ifdef debug_tips
			#pragma omp critical (cout)
						std::cerr << "MIN: " << min_max_walk_dist.first << " MAX: " << min_max_walk_dist.second << std::endl;
#endif
			std::vector<step_jaccard_t> target_jaccard_indices;
			// we were able to walk the given maximum walking distance in both directions, this greatly simplifies the algorithm
			if (min_max_walk_dist.first >= walking_dist && min_max_walk_dist.second >= walking_dist) {
				/// for the query:
				// walk the given walking_dist first to the left (we might not be able to do the full walk)
				// then walk the given walking_dist to the right
				// in order to collect all the visited nodes and how often they were visited in a ska::flat_hash_map<uint64_t, uint64_t>
				// the key is the node identifier, the value is the number of times this node was visited
				ska::flat_hash_map<nid_t, uint64_t> query_set = collect_nodes_in_walking_dist(graph, walking_dist, walking_dist, cur_step);
				// TODO we can precollect all sets so here we don't have to do the walks anymore
				for (step_handle_t& target_step : target_step_handles) {
					// we count which node identifier we saw how many times
					ska::flat_hash_map<nid_t, uint64_t> target_set = collect_nodes_in_walking_dist(graph, walking_dist, walking_dist, target_step);
					ska::flat_hash_map<nid_t, uint64_t> union_set;
					for (auto& query_item : query_set) {
						union_set[query_item.first] = query_item.second;
					}
#ifdef debug_tips
					std::cerr << "UNION SIZE BEFORE: " << union_set.size() << std::endl;
#endif
					add_target_set_to_union_set(union_set, target_set);
#ifdef debug_tips
					std::cerr << "UNION SIZE AFTER: " << union_set.size() << std::endl;


								for (auto u : union_set) {
									std::cerr << "node_id: " << u.first << " node_count: " << u.second << std::endl;
								}
#endif
					ska::flat_hash_map<nid_t, uint64_t> intersection_set = intersect_target_query_sets(union_set, target_set, query_set);
#ifdef debug_tips
					std::cerr << "INTERSECT SIZE: " << intersection_set.size() << std::endl;

								for (auto i : intersection_set) {
									std::cerr << "node_id: " << i.first << " node_count: " << i.second << std::endl;
								}
#endif
					double jaccard = jaccard_idx_from_intersect_union_sets(intersection_set, union_set, graph);
					target_jaccard_indices.push_back({target_step, jaccard});
#ifdef debug_tips
					#pragma omp critical (cout)
							std::cerr << "Jaccard index of query " << query_path_name << " and target " << target_path << ": " << jaccard << ". Direction: " << walk_from_front << std::endl;
#endif
				}
				// more complex algorithm
			} else {
				ska::flat_hash_map<nid_t, uint64_t> query_set_min_max = collect_nodes_in_walking_dist(
						graph,
						min_max_walk_dist.first,
						min_max_walk_dist.second,
						cur_step);
				ska::flat_hash_map<nid_t, uint64_t> query_set_max_min = collect_nodes_in_walking_dist(
						graph,
						min_max_walk_dist.second,
						min_max_walk_dist.first,
						cur_step);
				for (step_handle_t& target_step : target_step_handles) {
					ska::flat_hash_map<nid_t, uint64_t> target_set_min_max = collect_nodes_in_walking_dist(
							graph,
							min_max_walk_dist.first,
							min_max_walk_dist.second,
							target_step);
					ska::flat_hash_map<nid_t, uint64_t> target_set_max_min = collect_nodes_in_walking_dist(
							graph,
							min_max_walk_dist.second,
							min_max_walk_dist.first,
							target_step);

					/// [0] -> q_min_max vs. t_min_max
					/// [1] -> q_min_max vs. t_max_min
					/// [2] -> q_max_min vs. t_min_max
					/// [3] -> q_max_min vs. t_max_min
					/// will be 0.0 if a combinations is not possible
					std::vector<double> candidate_jaccards = {0.0, 0.0, 0.0, 0.0};
					if (query_set_min_max.size() > 0) {
						if (target_set_min_max.size() > 0) {
							candidate_jaccards[0] = get_jaccard_index(graph, query_set_min_max, target_set_min_max);
						}
						if (target_set_max_min.size() > 0) {
							candidate_jaccards[1] = get_jaccard_index(graph, query_set_min_max, target_set_max_min);
						}
					}
					if (query_set_max_min.size() > 0) {
						if (target_set_min_max.size() > 0) {
							candidate_jaccards[2] = get_jaccard_index(graph, query_set_max_min, target_set_min_max);
						}
						if (target_set_max_min.size() > 0) {
							candidate_jaccards[3] = get_jaccard_index(graph, query_set_max_min, target_set_max_min);
						}
					}
					target_jaccard_indices.push_back({target_step, *max_element(candidate_jaccards.begin(), candidate_jaccards.end())});
				}
			}
			std::sort(target_jaccard_indices.begin(), target_jaccard_indices.end(),
					  [&] (const step_jaccard_t & sjt_a,
						   const step_jaccard_t & sjt_b) {
						  return sjt_a.jaccard > sjt_b.jaccard;
					  });

			/// this ensures a deterministic selection of the best jaccard index and therefore step
			// collect all the step_jaccard_t with the same jaccard
			std::vector<step_jaccard_t> target_same_jaccard;
			bool first_pos = false;
			for (auto s_j_t : target_jaccard_indices) {
				if (first_pos) {
					if (s_j_t.jaccard == target_same_jaccard[target_same_jaccard.size() - 1].jaccard) {
						target_same_jaccard.push_back(s_j_t);
					} else {
						break;
					}
				} else {
					first_pos = true;
					target_same_jaccard.push_back(s_j_t);
				}
			}
			// sort by rank
			std::sort(target_same_jaccard.begin(), target_same_jaccard.end(),
					  [&] (const step_jaccard_t & sjt_a,
						   const step_jaccard_t & sjt_b) {
						  return as_integers(sjt_a.step)[1] < as_integers(sjt_b.step)[1];
					  });
			// take the one with array position arr_len/2; if arr_len%%2 !=0 then take the floor of the resulting value
			uint64_t final_jaccard_position_in_same = target_same_jaccard.size() / 2;
			step_jaccard_t final_jaccard = target_same_jaccard[final_jaccard_position_in_same];
			uint64_t final_jaccard_position;
			for (uint64_t i = 0; i < target_jaccard_indices.size(); i++) {
				step_jaccard_t s_j_t = target_jaccard_indices[i];
				if (s_j_t.jaccard == final_jaccard.jaccard
					&& as_integers(s_j_t.step)[0] == as_integers(final_jaccard.step)[0]
					&& as_integers(s_j_t.step)[1] == as_integers(final_jaccard.step)[1]) {
					final_jaccard_position = i;
					break;
				}
			}
			// use std::swap to switch the found array position with the first position
			std::swap(target_jaccard_indices[0], target_jaccard_indices[final_jaccard_position]);
			return target_jaccard_indices;
		}

		ska::flat_hash_map<nid_t , uint64_t> collect_nodes_in_walking_dist(const graph_t& graph,
																		   const uint64_t& walking_dist_prev,
																		   const uint64_t& walking_dist_next,
																		   const step_handle_t& start_step,
																		   const bool& walked_walking_dist) {
			ska::flat_hash_map<nid_t, uint64_t> node_count_set;
			/// first walk to previous steps up to the walking_dist
			uint64_t dist_walked = 0;
			uint64_t total_dist_walked = 0;
			// where does our step come from
			handle_t cur_h = graph.get_handle_of_step(start_step);
			step_handle_t cur_step = start_step;
			nid_t cur_id = graph.get_id(cur_h);
			while (graph.has_previous_step(cur_step) && (dist_walked < walking_dist_prev)) {
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
			total_dist_walked += dist_walked;
			dist_walked = 0;
			/// walking the next steps up to the walking_dist
			cur_step = start_step;
			while (graph.has_next_step(cur_step) && (dist_walked < walking_dist_next)) {
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
			total_dist_walked += dist_walked;
			// where does our step come from, we add id regardless of walking distance and orientation
			if (node_count_set.count(cur_id) == 0) {
				node_count_set[cur_id] = 1;
			} else {
				node_count_set[cur_id] = node_count_set[cur_id] + 1;
			}
			if ((total_dist_walked < (walking_dist_prev + walking_dist_next))) {
				return ska::flat_hash_map<nid_t, uint64_t>();
			}
			return node_count_set;
		}

		ska::flat_hash_map<nid_t , uint64_t> collect_nodes_in_walking_dist_from_map(const graph_t& graph,
																					const uint64_t& walking_dist_prev,
																					const uint64_t& walking_dist_next,
																					const step_handle_t& start_step,
																					ska::flat_hash_map<step_handle_t, std::pair<std::vector<nid_t>, std::vector<nid_t>>>& steps_nodes_prev_next_map) {
			ska::flat_hash_map<nid_t, uint64_t> node_count_set;
			/// first walk to previous steps up to the walking_dist
			uint64_t dist_walked = 0;
			uint64_t total_dist_walked = 0;
			// where does our step come from
			handle_t cur_h = graph.get_handle_of_step(start_step);
			nid_t cur_id = graph.get_id(cur_h);
			uint64_t idx = 0;
			vector<nid_t> prev_v = steps_nodes_prev_next_map[start_step].first;
			uint64_t prev_v_size = prev_v.size();
			while (idx < prev_v_size && (dist_walked < walking_dist_prev)) {
				nid_t prev_id = prev_v[idx];
				handle_t prev_h = graph.get_handle(prev_id);
				//
				if (node_count_set.count(prev_id) == 0) {
					node_count_set[prev_id] = 1;
				} else {
					node_count_set[prev_id] = node_count_set[prev_id] + 1;
				}
				dist_walked += graph.get_length(prev_h);
				idx++;
			}
			total_dist_walked += dist_walked;
			dist_walked = 0;
			/// walking the next steps up to the walking_dist
			idx = 0;
			vector<nid_t> next_v = steps_nodes_prev_next_map[start_step].second;
			uint64_t next_v_size = next_v.size();
			while (idx < next_v_size && dist_walked < walking_dist_next) {
				nid_t next_id = next_v[idx];
				handle_t next_h = graph.get_handle(next_id);
				if (node_count_set.count(next_id) == 0) {
					node_count_set[next_id] = 1;
				} else {
					node_count_set[next_id] = node_count_set[next_id] + 1;
				}
				dist_walked += graph.get_length(next_h);
				idx++;
			}
			total_dist_walked += dist_walked;
			// where does our step come from, we add id regardless of walking distance and orientation
			if (node_count_set.count(cur_id) == 0) {
				node_count_set[cur_id] = 1;
			} else {
				node_count_set[cur_id] = node_count_set[cur_id] + 1;
			}
			if ((total_dist_walked < (walking_dist_prev + walking_dist_next))) {
				return ska::flat_hash_map<nid_t, uint64_t>();
			}
			return node_count_set;
		}

		void add_target_set_to_union_set(ska::flat_hash_map<nid_t , uint64_t>& union_set,
								   const ska::flat_hash_map<nid_t , uint64_t>& target_set) {
			for (auto& target_item : target_set) {
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
			for (auto& union_item : union_set) {
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

		double get_jaccard_index(const graph_t& graph, ska::flat_hash_map<nid_t , uint64_t>& query_set,
								 ska::flat_hash_map<nid_t , uint64_t>& target_set) {
			ska::flat_hash_map<nid_t, uint64_t> union_set;
			for (auto& query_item : query_set) {
				union_set[query_item.first] = query_item.second;
			}
			add_target_set_to_union_set(union_set, target_set);
			ska::flat_hash_map<nid_t, uint64_t> intersection_set = intersect_target_query_sets(union_set, target_set, query_set);
			double jaccard = jaccard_idx_from_intersect_union_sets(intersection_set, union_set, graph);
			return jaccard;
		}

		double jaccard_idx_from_intersect_union_sets(ska::flat_hash_map<nid_t , uint64_t>& intersection_set,
													 ska::flat_hash_map<nid_t , uint64_t>& union_set,
													 const graph_t& graph) {
			double jaccard_idx = 0.0;
			uint64_t intersect_seq_len = 0;
			uint64_t union_seq_len = 0;
			for (auto& intersect_item : intersection_set) {
				nid_t intersect_node_id = intersect_item.first;
				uint64_t intersect_node_count = intersect_item.second;
				handle_t h = graph.get_handle(intersect_node_id);
				intersect_seq_len += graph.get_length(h)*intersect_node_count;
			}
			for (auto& union_item : union_set) {
				nid_t union_node_id = union_item.first;
				uint64_t union_node_count = union_item.second;
				handle_t h = graph.get_handle(union_node_id);
				union_seq_len += graph.get_length(h)*union_node_count;
			}
// #define debug_tips
#ifdef debug_tips
			std::cerr << "intersect_seq_len: " << intersect_seq_len << std::endl;
			std::cerr << "union_seq_len: " << union_seq_len << std::endl;
#endif
			jaccard_idx = (double) intersect_seq_len / (double) union_seq_len;
			return jaccard_idx;
		}

		std::pair<uint64_t , uint64_t> find_min_max_walk_dist_from_query_targets(const graph_t& graph,
																				 const uint64_t& walking_dist,
																				 const step_handle_t& cur_step,
																				 const std::vector<step_handle_t>& target_step_handles) {
			std::pair<uint64_t, uint64_t> min_max_walk_dist = {walking_dist, walking_dist};
			std::vector<step_handle_t> query_target_step_handles(target_step_handles);
			query_target_step_handles.push_back(cur_step);
			for (auto& start_step : query_target_step_handles) {
				/// first walk to previous steps up to the walking_dist
				uint64_t dist_walked_prev = 0;
				handle_t cur_h = graph.get_handle_of_step(start_step);
				step_handle_t cur_step = start_step;
				nid_t cur_id = graph.get_id(cur_h);
				while (graph.has_previous_step(cur_step) && (dist_walked_prev < min_max_walk_dist.second)) {
					step_handle_t prev_step = graph.get_previous_step(cur_step);
					handle_t prev_h = graph.get_handle_of_step(prev_step);
					nid_t prev_id = graph.get_id(prev_h);
					dist_walked_prev += graph.get_length(prev_h);
					cur_step = prev_step;
				}
				uint64_t dist_walked_next = 0;

				/// walking the next steps up to the walking_dist
				cur_step = start_step;
				while (graph.has_next_step(cur_step) && dist_walked_next < min_max_walk_dist.second) {
					step_handle_t next_step = graph.get_next_step(cur_step);
					handle_t next_h = graph.get_handle_of_step(next_step);
					nid_t next_id = graph.get_id(next_h);
					dist_walked_next += graph.get_length(next_h);
					cur_step = next_step;
				}
				min_max_walk_dist.first = std::min(std::min(dist_walked_prev, dist_walked_next), min_max_walk_dist.first);
				min_max_walk_dist.second = std::min(std::max(dist_walked_prev, dist_walked_next), min_max_walk_dist.second);
			}

			return min_max_walk_dist;
		}
	}
}
