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
					   const uint64_t& walking_dist,
					   const bool& report_additional_jaccards) {

			const std::string target_path = graph.get_path_name(target_path_t);

			std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_meter;
			if (progress) {
				std::string end;
				if (walk_from_front) {
					end = "front";
				} else {
					end = "back";
				}
				std::string progress_message = "[odgi::tips::walk_tips] BED Progress for target path '"
											   + graph.get_path_name(target_path_t) + "' walking from the " + end + " of all query paths:";
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
						std::vector<step_handle_t> target_step_handles;
						graph.for_each_step_on_handle(
								cur_h,
								[&](const step_handle_t& s) {
									/// we can do these expensive iterations here, because we only have to do it once for each walk
									if (graph.get_path_handle_of_step(s) == target_path_t) {
										// collect only the steps for the given target
										target_step_handles.push_back(s);
									}
								});

						std::vector<step_jaccard_t> target_jaccard_indices = jaccard_indices_from_step_handles(graph,
																											   walking_dist,
																											   cur_step,
																											   target_step_handles);
						uint64_t i = 0;

						// report other jaccards as a csv list in the BED
						std::vector<double> additional_jaccards_to_report;
						uint64_t start_index = 0 + n_best_mappings;
						if (!report_additional_jaccards) {
							start_index = target_jaccard_indices.size();
						}
						/// do we even have indices left for reporting?
						if (!(start_index >= target_jaccard_indices.size())) {
							for (uint64_t n = start_index; n < target_jaccard_indices.size(); n++) {
								additional_jaccards_to_report.push_back(target_jaccard_indices[n].jaccard);
							}
						}
						/// only report the Nth final steps
						for (auto& target_jaccard_index : target_jaccard_indices) {
							if (i == n_best_mappings) {
								break;
							}
							step_handle_t final_target_step = target_jaccard_index.step;
							double final_target_jaccard = target_jaccard_index.jaccard;

							uint64_t target_min_pos = step_index.get_position(final_target_step, graph); // 0-based starting position in BED
							uint64_t target_max_pos = target_min_pos + graph.get_length(graph.get_handle_of_step(final_target_step)); // 1-based ending position in BED

							/// add BED record to queue of the BED writer
							bed_writer_thread.append(target_path, target_min_pos, target_max_pos,
													 query_path_name, step_index.get_position(cur_step, graph),
													 final_target_jaccard, walk_from_front, additional_jaccards_to_report);
							i++;
						}
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
