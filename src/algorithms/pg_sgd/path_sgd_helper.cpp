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

		const void prepare_weak_connected_components_map(graph_t& graph,
											   std::vector<uint64_t>& weak_components_map) {
#ifdef debug_components
			std::cerr << "node count: " << graph.get_node_count() << std::endl;
#endif
			// refine order by weakly connected components

			// prepare weakly connected components
			std::vector<ska::flat_hash_set<handlegraph::nid_t>> weak_components = algorithms::weakly_connected_components(
					&graph);
#ifdef debug_components
			std::cerr << "components count: " << weak_components.size() << std::endl;
#endif
			std::vector<std::pair<double, uint64_t>> weak_component_order;
			for (int i = 0; i < weak_components.size(); i++) {
				auto &weak_component = weak_components[i];
				uint64_t id_sum = 0;
				for (auto node_id : weak_component) {
					id_sum += node_id;
				}
				double avg_id = id_sum / (double) weak_component.size();
				weak_component_order.push_back(std::make_pair(avg_id, i));
			}
			std::sort(weak_component_order.begin(), weak_component_order.end());
			std::vector<uint64_t> weak_component_id; // maps rank to "id" based on the original sorted order
			weak_component_id.resize(weak_component_order.size());
			uint64_t component_id = 0;
			for (auto &component_order : weak_component_order) {
				weak_component_id[component_order.second] = component_id++;
			}
			weak_components_map.resize(graph.get_node_count());
			// reserve the space we need
			for (int i = 0; i < weak_components.size(); i++) {
				auto &weak_component = weak_components[i];
				// store for each node identifier to component start index
				for (auto node_id : weak_component) {
					weak_components_map[node_id - 1] = weak_component_id[i];
				}
#ifdef debug_components
				std::cerr << "weak_component.size(): " << weak_component.size() << std::endl;
                std::cerr << "component_index: " << i << std::endl;
#endif
			}
			weak_components_map.clear();
		}

		/// prepare the weak components map which we will need later so we separate all graph components in 1D
		const void generate_and_write_snapshot_graphs(graph_t& graph,
													  std::vector<uint64_t>& weak_components_map,
													  const std::string& snapshot_prefix,
													  std::vector<string>& snapshots) {
			for (int j = 0; j < snapshots.size(); j++) {
				std::string snapshot_file_name = snapshots[j];
				std::ifstream snapshot_instream(snapshot_file_name);
				std::vector<double> snapshot_layout;
				std::string line;
				while(std::getline(snapshot_instream, line)) {
					snapshot_layout.push_back(std::stod(line));
				}
				snapshot_instream.close();
				uint64_t i = 0;
				std::vector<algorithms::handle_layout_t> snapshot_handle_layout;
				graph.for_each_handle(
						[&i, &snapshot_layout, &weak_components_map, &snapshot_handle_layout](
								const handle_t &handle) {
							snapshot_handle_layout.push_back(
									{
											weak_components_map[number_bool_packing::unpack_number(handle)],
											snapshot_layout[i++],
											handle
									});
						});
				// sort the graph layout by component, then pos, then handle rank
				std::sort(snapshot_handle_layout.begin(), snapshot_handle_layout.end(),
						  [&](const algorithms::handle_layout_t &a,
							  const algorithms::handle_layout_t &b) {
							  return a.weak_component < b.weak_component
									 || (a.weak_component == b.weak_component
										 && a.pos < b.pos
										 || (a.pos == b.pos
											 && as_integer(a.handle) < as_integer(b.handle)));
						  });
				std::vector<handle_t> order;
				order.reserve(graph.get_node_count());
				for (auto &layout_handle : snapshot_handle_layout) {
					order.push_back(layout_handle.handle);
				}
				std::cerr << "[odgi::path_linear_sgd] Applying order to graph of snapshot: " << std::to_string(j + 1)
						  << std::endl;
				std::string local_snapshot_prefix = snapshot_prefix + std::to_string(j + 1);
				auto* graph_copy = new odgi::graph_t();
				utils::graph_deep_copy(graph, graph_copy);
				graph_copy->apply_ordering(order, true);
				ofstream f(local_snapshot_prefix);
				std::cerr << "[odgi::path_linear_sgd] Writing snapshot: " << std::to_string(j + 1) << std::endl;
				graph_copy->serialize(f);
				f.close();
			}
		}

		const void from_layout_to_node_order(graph_t& graph,
											 std::vector<uint64_t>& weak_components_map,
											 std::vector<handle_t>& order,
											 std::vector<double>& layout) {
			std::vector<algorithms::handle_layout_t> handle_layout;
			uint64_t i = 0;
			graph.for_each_handle(
					[&i, &layout, &weak_components_map, &handle_layout](const handle_t &handle) {
						handle_layout.push_back(
								{
										weak_components_map[number_bool_packing::unpack_number(handle)],
										layout[i++],
										handle
								});
					});
			// sort the graph layout by component, then pos, then handle rank
			std::sort(handle_layout.begin(), handle_layout.end(),
					  [&](const algorithms::handle_layout_t &a,
						  const algorithms::handle_layout_t &b) {
						  return a.weak_component < b.weak_component
								 || (a.weak_component == b.weak_component
									 && a.pos < b.pos
									 || (a.pos == b.pos
										 && as_integer(a.handle) < as_integer(b.handle)));
					  });
			order.reserve(graph.get_node_count());
			for (auto &layout_handle : handle_layout) {
				order.push_back(layout_handle.handle);
			}
		}
	}
}
