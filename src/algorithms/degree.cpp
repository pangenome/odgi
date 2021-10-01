#include "degree.hpp"

namespace odgi {
	namespace algorithms {

		void for_each_path_range_degree(const PathHandleGraph& graph,
									   const std::vector<path_range_t>& _path_ranges,
									   const std::vector<bool>& paths_to_consider,
									   const std::function<void(const path_range_t&, const double&)>& func) {
			const uint64_t shift = graph.min_node_id();
			if (graph.max_node_id() - shift >= graph.get_node_count()){
				std::cerr << "[degree::for_each_path_range_degree] error: the node IDs are not compacted. Please run 'odgi sort' using -O, --optimize to optimize the graph." << std::endl;
				exit(1);
			}
			// are we subsetting?
			const bool subset_paths = !paths_to_consider.empty();
			// sort path ranges
			std::vector<const path_range_t*> path_ranges;
			for (auto& p : _path_ranges) {
				path_ranges.push_back(&p);
			}
			std::sort(path_ranges.begin(),
					  path_ranges.end(),
					  [&](const path_range_t* a,
						  const path_range_t* b) {
						  return std::tie(as_integer(a->begin.path), a->begin.offset, a->end.offset, a->begin.is_rev)
								 < std::tie(as_integer(b->begin.path), b->begin.offset, b->end.offset, b->begin.is_rev);
					  });
			// order our path ranges
			std::vector<std::pair<std::vector<const path_range_t*>::const_iterator,
					std::vector<const path_range_t*>::const_iterator>> paths_todo;
			auto p = path_ranges.begin();
			while (p != path_ranges.end()) {
				auto n = p; ++n;
				while (n != path_ranges.end()
					   && (*n)->begin.path == (*p)->begin.path) {
					++n;
				}
				paths_todo.push_back(std::make_pair(p, n));
				p = n;
			}
			// precompute degrees for all handles in parallel
			std::vector<uint64_t> degree(graph.get_node_count() + 1);
			graph.for_each_handle(
					[&](const handle_t& h) {
						auto id = graph.get_id(h);
						auto& d = degree[id - shift];
						if (subset_paths) {
							bool consider = false;
							graph.for_each_step_on_handle(
									h,
									[&](const step_handle_t &s) {
										if (paths_to_consider[as_integer(graph.get_path_handle_of_step(s))]) {
											consider = true;
										}
									});
							if (consider) {
								d = graph.get_degree(h, false) + graph.get_degree(h, true);
							}
						} else {
							d = graph.get_degree(h, false) + graph.get_degree(h, true);
						}
					}, true);
			// dip into the path ranges, for each
#pragma omp parallel for schedule(dynamic,1)
			for (auto& p : paths_todo) {
				auto& path = (*p.first)->begin.path;
				std::vector<std::pair<const path_range_t*, uint64_t>> active_ranges;
				auto range_itr = p.first;
				auto ranges_end = p.second;
				uint64_t offset = 0;
				graph.for_each_step_in_path(
						path,
						[&](const step_handle_t& step) {
							// understand our context
							handle_t handle = graph.get_handle_of_step(step);
							auto node_length = graph.get_length(handle);
							auto node_end = offset + node_length;
							// remove non-active ranges and call callbacks
							active_ranges.erase(
									std::remove_if(active_ranges.begin(),
												   active_ranges.end(),
												   [&](auto& r) {
													   if (r.first->end.offset <= offset) {
														   func(*r.first,
																(double)r.second /
																(double)(r.first->end.offset
																		 - r.first->begin.offset));
														   return true;
													   } else {
														   return false;
													   }
												   }),
									active_ranges.end());
							// find ranges that start in this node, and add them to active ranges
							while (range_itr != ranges_end
								   && (*range_itr)->begin.offset < node_end
								   && (*range_itr)->begin.offset >= offset) {
								active_ranges.push_back(std::make_pair(*range_itr++, 0));
							}
							// compute degree on node
							auto& d = degree[graph.get_id(handle) - shift];
							// add coverage to active ranges
							for (auto& r : active_ranges) {
								auto l = node_length;
								if (r.first->begin.offset > offset) {
									l -= r.first->begin.offset - offset;
								}
								if (r.first->end.offset < node_end) {
									l -= node_end - r.first->end.offset;
								}
								r.second += (d * l);
							}
							// to node end
							offset = node_end;
						});
				// clear out any ranges that are still active at the end
				for (auto& r : active_ranges) {
					func(*r.first,
						 (double)r.second /
						 (double)(r.first->end.offset
								  - r.first->begin.offset));
				}
			}
		}

	}
}
