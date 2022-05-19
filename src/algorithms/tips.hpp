#pragma once

#include "algorithms/stepindex.hpp"
#include "algorithms/tips_bed_writer_thread.hpp"
#include "algorithms/path_jaccard.hpp"
#include "odgi.hpp"
#include <omp.h>
#include "hash_map.hpp"

/**
 * \file tips.hpp
 *
 * Defines algorithms for identifying break point positions relative to given query (reference) path(s) of all the tips in the graph or of tips of given path(s).
 */

namespace odgi {
	namespace algorithms {

		using namespace handlegraph;

		/// Iterate over the given paths. We walk from a chosen end (path_start or path_back) until we hit a node of the given query (reference) path(s).
		/// We record the hit as a BED output.
		/// #chrom #start #end #path_name #path_pos #jaccard #from_front #add_jaccards
		/// #chrom: The query path name.
		/// #start: The 0-based start position of the query we hit in the node.
		/// #end: The 1-based end position of the query we hit in the node.
		/// #path: The name of the path we walked.
		/// #path_pos: The 0-based position of the path we walked when we hit the node of the query path.
		/// #jaccard: The jaccard index of the query and target path around the region of the step where the query hit the target.
		/// #walk_from_front: If 1 we walked from the head of the target path. Else we walked from the tail and it is 0.
		/// add_jaccards: The additional jaccards of candidate reference step(s). Comma-separated.
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
				 const bool& report_additional_jaccards);
	}
}
