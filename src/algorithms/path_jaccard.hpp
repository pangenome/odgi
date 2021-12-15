#pragma once

#include "algorithms/stepindex.hpp"
#include "algorithms/tips_bed_writer_thread.hpp"
#include <omp.h>
#include "hash_map.hpp"

/**
 * \file path_jaccard.hpp
 *
 * Defines algorithms for the calculation of the path jaccard concept, an alignment-free local sequence similarity threshold.
 */

namespace odgi {
	namespace algorithms {

		using namespace handlegraph;

		/// here we store a jaccard index together with the corresponding step
		struct step_jaccard_t {
			step_handle_t step;
			double jaccard = 0.0;
		};

		/// calculate all jaccard indices from a given target_step_handles and a current query step
		/// the MAJOR function!
		std::vector<step_jaccard_t> jaccard_indices_from_step_handles(const graph_t& graph,
																	  const uint64_t& walking_dist,
																	  const step_handle_t& cur_step,
																	  std::vector<step_handle_t>& target_step_handles);

		/// from the given start step we walk the given distance in nucleotides left and right following the steps in the given graph graph, collecting all nodes that we cross <key>
		/// we also record, how many times we visited a node <value>
		ska::flat_hash_map<nid_t , uint64_t> collect_nodes_in_walking_dist(const graph_t& graph,
																		const uint64_t& walking_dist_prev,
																		const uint64_t& walking_dist_next,
																		const step_handle_t& start_step,
																		const bool& walked_walking_dist = false);

		/// from the give start step we walk the given distance in nucleotides left and right following the given map, collecting all nodes that we cross <key>
		/// we also record, how many times we visited a node <value>
		ska::flat_hash_map<nid_t , uint64_t> collect_nodes_in_walking_dist_from_map(const graph_t& graph,
																		   const uint64_t& walking_dist_prev,
																		   const uint64_t& walking_dist_next,
																		   const step_handle_t& start_step);

		/// from a given target_set add the nodes into the union_set which might be not empty
		void add_target_set_to_union_set(ska::flat_hash_map<nid_t , uint64_t>& union_set,
								   const ska::flat_hash_map<nid_t , uint64_t>& target_set);

		/// from a given target_set and a given query_set generate the intersection_set
		/// we use the union_set as guidance
		ska::flat_hash_map<nid_t, uint64_t> intersect_target_query_sets(ska::flat_hash_map<nid_t , uint64_t>& union_set,
										 ska::flat_hash_map<nid_t , uint64_t>& target_set,
										 ska::flat_hash_map<nid_t , uint64_t>& query_set);

		/// calculate the jaccard index from a given graph, a query set and a target set
		/// invokes jaccard_idx_from_intersect_union_sets
		double get_jaccard_index(const graph_t& graph, ska::flat_hash_map<nid_t , uint64_t>& query_set,
								 ska::flat_hash_map<nid_t , uint64_t>& target_set);

		/// calculate the jaccard index from an intersection_set and a union_set
		/// 1. calculate the sequence lengths of both sets
		/// 2. intersection_set_seq_len / union_set_seq_len
		double jaccard_idx_from_intersect_union_sets(ska::flat_hash_map<nid_t , uint64_t>& intersection_set,
															   ska::flat_hash_map<nid_t , uint64_t>& union_set,
															   const graph_t& graph);

		/// given a vector of target step handles and a walking distance, we want to find out how much of the walking distance we can follow into each direction for each step
		/// we identify the set of the maximum walkable distance for both directions shared by all steps
		std::pair<uint64_t , uint64_t> find_min_max_walk_dist_from_query_targets(const graph_t& graph,
																		   const uint64_t& walking_dist,
																		   const step_handle_t& cur_step,
																		   const std::vector<step_handle_t>& target_step_handles);
	}
}
