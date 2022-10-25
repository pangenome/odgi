#pragma once

#include "odgi.hpp"
#include "xp.hpp"
#include "algorithms/stepindex.hpp"

namespace odgi {
	namespace algorithms {

		using namespace handlegraph;

		/// get the total number of steps in the graph
		const uint64_t get_sum_path_step_count(const std::vector<path_handle_t> &path_sgd_use_paths,
											   graph_t &graph);

		/// get the maximum step count of all paths in the graph
		const uint64_t get_max_path_step_count(const std::vector<path_handle_t> &path_sgd_use_paths,
											   graph_t &graph);

		/// get the maximum path length in nucleotides using the XP
		const uint64_t get_max_path_length_xp(const std::vector<path_handle_t> &path_sgd_use_paths,
											  const xp::XP &path_index);

		/// get the maximum path length in nucleotides using the SSI
		const uint64_t get_max_path_length_ssi(const std::vector<path_handle_t> &path_sgd_use_paths,
											   const algorithms::step_index_t &sampled_step_index);
	}
}