#include "tips.hpp"

namespace odgi {
	namespace algorithms {

		using namespace handlegraph;

		void walk_tips(const graph_t& graph,
					   const std::vector<path_handle_t>& paths,
					   const path_handle_t& query_path_t,
					   const std::vector<bool>& query_handles,
					   algorithms::step_index_t& step_index,
					   const uint64_t& num_threads,
					   const std::function<step_handle_t(const path_handle_t&)>& get_path_end,
					   const std::function<step_handle_t(const step_handle_t&)>& get_step,
					   algorithms::tips_bed_writer& bed_writer_thread,
					   const bool progress,
					   const bool walk_from_front) {

			std::string query_path = graph.get_path_name(query_path_t);

			std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_meter;
			if (progress) {
				std::string end;
				if (walk_from_front) {
					end = "front";
				} else {
					end = "back";
				}
				std::string progress_message = "[odgi::tips::walk_tips] BED Progress for query path '"
						+ graph.get_path_name(query_path_t) + "' walking from the " + end + " of all paths:";
				progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
						paths.size(), progress_message);
			}

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
			for (auto path : paths) {
				// prevent self tips
				if (path == query_path_t) {
					continue;
				}
				std::string target_path_name = graph.get_path_name(path);
				bool tip_reached_query = false;
				// collecting tips from the front
				step_handle_t cur_step = get_path_end(path);
				handle_t cur_h = graph.get_handle_of_step(cur_step);
				while (!tip_reached_query) {
					// did we already hit the given reference path?
					if (query_handles[number_bool_packing::unpack_number(cur_h)]) {
						std::vector<uint64_t> query_path_positions;
						graph.for_each_step_on_handle(
								cur_h,
								[&](const step_handle_t& s) {
									if (graph.get_path_handle_of_step(s) == query_path_t) {
//#pragma omp critical (cout)
//									std::cerr << "[odgi::tips] error: We reached a tip from the front!" << std::endl;

										// we collect all positions
										// later we will get the smallest, the highest and the median from these
										query_path_positions.push_back(step_index.get_position(s));
									}
								});
//#pragma omp critical (cout)
//					std::cerr << "query_path_positions.size(): " << query_path_positions.size() << std::endl;
						std::sort(query_path_positions.begin(),
								  query_path_positions.end(),
								  [&](const uint64_t & pos_a,
									  const uint64_t & pos_b) {
									  return pos_a < pos_b;
								  });
						double query_pos_median = utils::median_of_sorted_vec(query_path_positions);
						uint64_t query_min_pos = query_path_positions[0]; // 0-based starting position in BED
						uint64_t query_max_pos = query_path_positions[query_path_positions.size() - 1] + 1; // 1-based ending position in BED
//#pragma omp critical (cout)
						//std::cout << query_path << "\t" << query_min_pos << "\t" << query_max_pos << "\t"
						//		  << query_pos_median << "\t" << graph.get_path_name(path) << "\t" << step_index.get_position(cur_step) << std::endl;
						/// add BED record to queue of BED_writer_thread
						bed_writer_thread.append(query_path, query_min_pos, query_max_pos,
							   query_pos_median, target_path_name, step_index.get_position(cur_step));
						tip_reached_query = true;
					}
					if (graph.has_next_step(cur_step)) {
						cur_step = get_step(cur_step);
						cur_h = graph.get_handle_of_step(cur_step);
					} else {
						// did we iterate over all steps and we did not hit the reference?
						tip_reached_query = true;
						// TODO emit a warning here on stderr
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
