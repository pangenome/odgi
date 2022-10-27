#pragma once

#include <vector>
#include <string>
#include <cassert>
#include <algorithm>
#include <random>
#include <set>
#include <thread>
#include <atomic>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/handle_graph.hpp>
#include "../stepindex.hpp"
#include "../xp.hpp" // we still ned the xp for the tmp_file, maybe move this to utils?
#include "../sgd_term.hpp"
#include "IITree.h"
#include <iomanip>
#include <string>
#include "../weakly_connected_components.hpp"
#include <sdsl/bit_vectors.hpp>
#include "dirty_zipfian_int_distribution.h"
#include "XoshiroCpp.hpp"
#include "../progress.hpp"
#include "utils.hpp"
#include "path_sgd_helper.hpp"

#include <fstream>

namespace odgi {
	namespace algorithms {

		using namespace handlegraph;

/// use SGD driven, by path guided, and partly zipfian distribution sampled pairwise distances to obtain a 1D linear layout of the graph that respects its topology
// FIXME make graph const again
		std::vector<double> path_linear_sgd_ssi(graph_t &graph,
											const algorithms::step_index_t &sampled_step_index,
											const std::vector<path_handle_t>& path_sgd_use_paths,
											const uint64_t &iter_max,
											const uint64_t &iter_with_max_learning_rate,
											const uint64_t &min_term_updates,
											const double &delta,
											const double &eps,
											const double &eta_max,
											const double &theta,
											const uint64_t &space,
											const uint64_t &space_max,
											const uint64_t &space_quantization_step,
											const double &cooling_start,
											const uint64_t &nthreads,
											const bool &progress,
											const bool &snapshot,
											std::vector<std::string> &snapshots,
											const uint64_t &sum_path_step_count);

/// our learning schedule
		std::vector<double> path_linear_sgd_schedule_ssi(const double &w_min,
													 const double &w_max,
													 const uint64_t &iter_max,
													 const uint64_t &iter_with_max_learning_rate,
													 const double &eps);
// FIMXME make grpah const again
		std::vector<handle_t> path_linear_sgd_order_ssi(graph_t &graph,
													const algorithms::step_index_t &sampled_step_index,
													const std::vector<path_handle_t>& path_sgd_use_paths,
													const uint64_t &iter_max,
													const uint64_t &iter_with_max_learning_rate,
													const uint64_t &min_term_updates,
													const double &delta,
													const double &eps,
													const double &eta_max,
													const double &theta,
													const uint64_t &space,
													const uint64_t &space_max,
													const uint64_t &space_quantization_step,
													const double &cooling_start,
													const uint64_t &nthreads,
													const bool &progress,
													const std::string &seed,
													const bool &snapshot,
													const std::string &snapshot_prefix,
													const bool &target_sorting,
													std::vector<bool>& target_nodes,
													const uint64_t &sum_path_step_count);

	}

}