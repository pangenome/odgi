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
	}
}
