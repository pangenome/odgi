#pragma once

#include "odgi.hpp"
#include "xp.hpp"
#include "algorithms/stepindex.hpp"

namespace odgi {
	namespace algorithms {

		using namespace handlegraph;

		struct handle_layout_t {
			uint64_t weak_component = 0;
			double pos = 0;
			handle_t handle = as_handle(0);
		};

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

		/// sort the graph's nodes by given target paths, we need this for the reference based sorting
		const void sort_graph_by_target_paths(graph_t& graph, const std::vector<path_handle_t> &target_paths,
											  std::vector<bool>& is_ref,
											  const bool &progress);

		/// calculate the schedule
		std::vector<double> path_linear_sgd_schedule(const double &w_min,
														 const double &w_max,
														 const uint64_t &iter_max,
														 const uint64_t &iter_with_max_learning_rate,
														 const double &eps);
	}
}
