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
						std::vector<uint64_t> target_path_start_positions;
						std::vector<uint64_t> target_path_end_positions;
						graph.for_each_step_on_handle(
								cur_h,
								[&](const step_handle_t& s) {
									/// we can do these expensive iterations here, because we only have to do it once for each walk
									if (graph.get_path_handle_of_step(s) == target_path_t) {
										// TODO collect only the steps for the given target
										// we collect all positions
										// later we will get the smallest, the highest and the median from these
										uint64_t pos = step_index.get_position(s);
										target_path_start_positions.push_back(pos);
										target_path_end_positions.push_back(pos + graph.get_length(cur_h) - 1); // 0-based
									}
								});
						/// for the query:
						// walk the given walking_dist first to the left (we might not be able to do the full walk)
						// then walk the given walking_dist to the right
						// in order to collect all the visited nodes and how often they were visited in a ska::flat_hash_map<uint64_t, uint64_t>
						// the key is the node identifier, the value is the number of times this node was visited
						// TODO pack the following into a function returning the query set, so we can use it for the query and all target's steps that we found
#pragma omp critical (cout)
						std::cerr << "query name: " << query_path_name << std::endl;
						ska::flat_hash_map<uint64_t, uint64_t> query_set;
						/// first walk to previous steps up to the walking_dist
						uint64_t dist_walked = 0;
						// where does our step come from
						const bool is_rev = number_bool_packing::unpack_bit(cur_h);
						step_handle_t query_cur_step = cur_step;
						nid_t cur_id = graph.get_id(cur_h);
						if (!is_rev) {
							dist_walked += graph.get_length(cur_h);
							query_set[cur_id] = 1;
						}
						while (graph.has_previous_step(query_cur_step) && (dist_walked <= walking_dist)) {
							step_handle_t prev_step = graph.get_previous_step(query_cur_step);
							handle_t prev_h = graph.get_handle_of_step(prev_step);
							nid_t prev_id = graph.get_id(prev_h);
							//
							if (query_set.count(prev_id) == 0) {
								query_set[prev_id] = 1;
							} else {
								query_set[prev_id] = query_set[prev_id] + 1;
							}
							dist_walked += graph.get_length(prev_h);
							query_cur_step = prev_step;
						}
						uint64_t dist_to_walk = walking_dist - dist_walked;
						dist_to_walk = (dist_to_walk > 0 ? dist_to_walk : 0);
						dist_walked = 0 - dist_to_walk;
						// where does our step come from
						if (is_rev) {
							dist_walked += graph.get_length(cur_h);
							if (query_set.count(cur_id) == 0) {
								query_set[cur_id] = 1;
							} else {
								query_set[cur_id] = query_set[cur_id] + 1;
							}
						}
						/// walking the next steps up to the walking_dist
						query_cur_step = cur_step;
						while (graph.has_next_step(query_cur_step) && dist_walked <= walking_dist) {
								step_handle_t next_step = graph.get_next_step(query_cur_step);
								handle_t next_h = graph.get_handle_of_step(next_step);
								nid_t next_id = graph.get_id(next_h);
								if (query_set.count(next_id) == 0) {
									query_set[next_id] = 1;
								} else {
									query_set[next_id] = query_set[next_id] + 1;
								}
								dist_walked += graph.get_length(next_h);
								query_cur_step = next_step;
						}
						for (auto pair : query_set) {
#pragma omp critical (cout)
							std::cerr << "first: " << pair.first << " second: " << pair.second << std::endl;
						}

						// TODO wrap the above into a function!!!
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

	}
}
