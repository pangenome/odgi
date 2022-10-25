#include "path_sgd_helper.hpp"

namespace odgi {
	namespace algorithms {

		using namespace handlegraph;

		const uint64_t get_sum_path_step_count(const std::vector<path_handle_t> &path_sgd_use_paths, graph_t &graph) {
			uint64_t sum_path_step_count = 0;
			for (auto& path : path_sgd_use_paths) {
				sum_path_step_count += graph.get_step_count(path);
			}
			return sum_path_step_count;
		};

		const uint64_t get_max_path_step_count(const std::vector<path_handle_t> &path_sgd_use_paths, graph_t &graph) {
			uint64_t max_path_step_count = 0;
			for (auto& path : path_sgd_use_paths) {
				max_path_step_count = std::max(max_path_step_count, graph.get_step_count(path));
			}
			return max_path_step_count;
		}

		const uint64_t get_max_path_length_xp(const std::vector<path_handle_t> &path_sgd_use_paths, const xp::XP &path_index) {
			uint64_t max_path_length = std::numeric_limits<uint64_t>::min();
			for (auto &path : path_sgd_use_paths) {
				max_path_length = std::max(max_path_length, path_index.get_path_length(path));
			}
			return max_path_length;
		}

		const uint64_t get_max_path_length_ssi(const std::vector<path_handle_t> &path_sgd_use_paths, const algorithms::step_index_t &sampled_step_index) {
			uint64_t max_path_length = std::numeric_limits<uint64_t>::min();
			for (auto &path : path_sgd_use_paths) {
				max_path_length = std::max(max_path_length, sampled_step_index.get_path_len(path));
			}
			return max_path_length;
		}

		const void sort_graph_by_target_paths(graph_t& graph, const std::vector<path_handle_t> &target_paths,
											  std::vector<bool>& is_ref,
											  const bool &progress) {
			std::vector<handle_t> target_order;
			std::fill_n(std::back_inserter(is_ref), graph.get_node_count(), false);
			std::unique_ptr <odgi::algorithms::progress_meter::ProgressMeter> target_paths_progress;
			if (progress) {
				std::string banner = "[odgi::sort] preparing target path vectors:";
				target_paths_progress = std::make_unique<odgi::algorithms::progress_meter::ProgressMeter>(target_paths.size(), banner);
			}
			for (handlegraph::path_handle_t target_path: target_paths) {
				graph.for_each_step_in_path(
						target_path,
						[&](const step_handle_t &step) {
							handle_t handle = graph.get_handle_of_step(step);
							uint64_t i = graph.get_id(handle) - 1;
							if (!is_ref[i]) {
								is_ref[i] = true;
								target_order.push_back(handle);
							}
						});
				if (progress) {
					target_paths_progress->increment(1);
				}
			}
			if (progress)  {
				target_paths_progress->finish();
			}
			uint64_t ref_nodes = 0;
			for (uint64_t i = 0; i < is_ref.size(); i++) {
				bool ref = is_ref[i];
				if (!ref) {
					target_order.push_back(graph.get_handle(i + 1));
					ref_nodes++;
				}
			}
			graph.apply_ordering(target_order, true);

			// refill is_ref with start->ref_nodes: 1 and ref_nodes->end: 0
			std::fill_n(is_ref.begin(), ref_nodes, true);
			std::fill(is_ref.begin() + ref_nodes, is_ref.end(), false);
		}
	}
}
