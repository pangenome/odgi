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

		std::vector<double> path_linear_sgd_schedule(const double &w_min,
														 const double &w_max,
														 const uint64_t &iter_max,
														 const uint64_t &iter_with_max_learning_rate,
														 const double &eps) {
#ifdef debug_schedule
			std::cerr << "w_min: " << w_min << std::endl;
            std::cerr << "w_max: " << w_max << std::endl;
            std::cerr << "iter_max: " << iter_max << std::endl;
            std::cerr << "eps: " << eps << std::endl;
#endif
			double eta_max = 1.0 / w_min;
			double eta_min = eps / w_max;
			double lambda = log(eta_max / eta_min) / ((double) iter_max - 1);
#ifdef debug_schedule
			std::cerr << "eta_max: " << eta_max << std::endl;
            std::cerr << "eta_min: " << eta_min << std::endl;
            std::cerr << "lambda: " << lambda << std::endl;
#endif
			// initialize step sizes
			std::vector<double> etas;
			etas.reserve(iter_max + 1);
#ifdef debug_schedule
			std::cerr << "etas: ";
#endif
			for (int64_t t = 0; t <= iter_max; t++) {
				etas.push_back(eta_max * exp(-lambda * (abs(t - (int64_t) iter_with_max_learning_rate))));
#ifdef debug_schedule
				std::cerr << etas.back() << ", ";
#endif
			}
#ifdef debug_schedule
			std::cerr << std::endl;
#endif
			return etas;
		}
	}
}
