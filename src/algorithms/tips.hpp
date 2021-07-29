#pragma once

#include "utils.hpp"
#include "algorithms/stepindex.hpp"
#include "algorithms/tips_bed_writer_thread.hpp"
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
		/// #chrom #start #end #median_range #path_name #path_pos
		/// #chrom: The query path name.
		/// #start: The 0-based start position of the query we hit in the node.
		/// #end: The 1-based end position of the query we hit in the node.
		/// #median_range: The 0-based median of the whole query path range of the node we hit. It is possible that a node contains several steps, so we want to mirror that here.
		/// #path: The name of the path we walked.
		/// #path_pos: The 0-based position of the path we walked when we hit the node of the query path.
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

		/// from the given start_step we walke the given distance in nucleotides left and right, collecting all nodes that we cross <key>
		/// we also record, how many times we visited a node <value>
		ska::flat_hash_map<nid_t , uint64_t> collect_nodes_in_walking_dist(const graph_t& graph,
																		const uint64_t& walking_dist,
																		const step_handle_t& start_step);

		/// from a given target_set add the nodes into the union_set which might be not empty
		void add_target_set_to_union_set(ska::flat_hash_map<nid_t , uint64_t>& union_set,
								   const ska::flat_hash_map<nid_t , uint64_t>& target_set);

		/// from a given target_set and a given query_set generate the intersection_set
		/// we use the union_set as guidance
		ska::flat_hash_map<nid_t, uint64_t> intersect_target_query_sets(ska::flat_hash_map<nid_t , uint64_t>& union_set,
										 ska::flat_hash_map<nid_t , uint64_t>& target_set,
										 ska::flat_hash_map<nid_t , uint64_t>& query_set);

		/// calculate the jaccard index from an intersection_set and a union_set
		/// 1. calculate the sequence lengths of both sets
		/// 2. intersection_set_seq_len / union_set_seq_len
		double jaccard_idx_from_intersect_union_sets(ska::flat_hash_map<nid_t , uint64_t>& intersection_set,
															   ska::flat_hash_map<nid_t , uint64_t>& union_set,
															   const graph_t& graph);
		struct step_jaccard_t {
			step_handle_t step;
			double jaccard = 0.0;
		};
	}
}
