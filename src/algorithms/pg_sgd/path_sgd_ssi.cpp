#include "path_sgd_ssi.hpp"

//#define debug_path_sgd
// #define eval_path_sgd
// #define debug_schedule
// # define debug_sample_from_nodes
namespace odgi {
	namespace algorithms {

		std::vector<double> path_linear_sgd_ssi(graph_t &graph,
											const algorithms::step_index_t &sampled_step_index,
											const std::vector<path_handle_t> &path_sgd_use_paths,
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
											const bool &target_sorting,
											std::vector<bool>& target_nodes,
											const uint64_t &sum_path_step_count) {
#ifdef debug_path_sgd
			std::cerr << "iter_max: " << iter_max << std::endl;
            std::cerr << "min_term_updates: " << min_term_updates << std::endl;
            std::cerr << "delta: " << delta << std::endl;
            std::cerr << "eps: " << eps << std::endl;
            std::cerr << "theta: " << theta << std::endl;
            std::cerr << "space: " << space << std::endl;
            std::cerr << "space_max: " << space_max << std::endl;
            std::cerr << "space_quantization_step: " << space_quantization_step << std::endl;
            std::cerr << "cooling_start: " << cooling_start << std::endl;
#endif

			uint64_t first_cooling_iteration = std::floor(cooling_start * (double)iter_max);
			//std::cerr << "first cooling iteration " << first_cooling_iteration << std::endl;

			/*
			// FIXME this is for testing only
			sdsl::bit_vector s_bv;
			sdsl::util::assign(s_bv, sdsl::bit_vector(10));
			s_bv[0] = 1;
			s_bv[3] = 1;
			s_bv[5] = 1;
			s_bv[8] = 1;
			for (auto bitty : s_bv) {
				std::cerr << bitty << std::endl;
			}
			sdsl::rank_support_v<1> s_bv_rank;
			sdsl::bit_vector::select_1_type s_bv_select;
			sdsl::util::assign(s_bv_rank, sdsl::rank_support_v<1>(&s_bv));
			sdsl::util::assign(s_bv_select, sdsl::bit_vector::select_1_type(&s_bv));
			std::cerr << "assignment complete" << std::endl;
			uint64_t k = 6;
			uint64_t k_rank = s_bv_rank(k);
			std::cerr << "rank at 6: " << k_rank << std::endl;
			uint64_t k_rank_to_select = s_bv_select(k_rank);
			std::cerr << "select at rank " << k_rank << ": " << k_rank_to_select << std::endl;
			uint64_t step_on_node = k - k_rank_to_select - 1; // 0-based!
			std::cerr << "step on node 0-based: " << step_on_node << std::endl;
			// this case will be skipped anyhow
			k = 5;
			k_rank = s_bv_rank(k);
			std::cerr << "rank at 6: " << k_rank << std::endl;
			k_rank_to_select = s_bv_select(k_rank);
			std::cerr << "select at rank " << k_rank << ": " << k_rank_to_select << std::endl;
			step_on_node = k - k_rank_to_select - 1; // 0-based!
			std::cerr << "step on node 0-based: " << step_on_node << std::endl;
			*/

			/// initialize the bitvector for random sampling
			/// each time we see a node we add 1, and each time we see a step we keep it 0
			sdsl::bit_vector ns_bv; // node-step bitvector
			sdsl::util::assign(ns_bv, graph.get_node_count() + sum_path_step_count);
			/// to be super save we also record the orientation of each step
			uint64_t len = 0;
			uint64_t num_nodes = graph.get_node_count();
			// our positions in 1D
			std::vector<std::atomic<double>> X(num_nodes);
			std::vector<bool> handle_orientation;
			handle_orientation.reserve(sum_path_step_count);
			uint64_t l = 0;
			uint64_t m = 0;

			// add progress bar here?
			graph.for_each_handle(
					[&](const handle_t &h) {
						ns_bv[l] = 1;
						// nb: we assume that the graph provides a compact handle set
						X[number_bool_packing::unpack_number(h)].store(len);
						len += graph.get_length(h);
						graph.for_each_step_on_handle(
								h,
								[&](const step_handle_t& step) {
									l++;
									handle_orientation[m] = graph.get_is_reverse(graph.get_handle_of_step(step));
									m++;
								});
						l++;
					});
			sdsl::rank_support_v<1> ns_bv_rank;
			sdsl::bit_vector::select_1_type ns_bv_select;
			sdsl::util::assign(ns_bv_rank, sdsl::rank_support_v<1>(&ns_bv));
			sdsl::util::assign(ns_bv_select, sdsl::bit_vector::select_1_type(&ns_bv));

			/// FIXME this is just so that it runs
			xp::XP path_index;
			path_index.from_handle_graph(graph, nthreads);

			uint64_t total_term_updates = iter_max * min_term_updates;
			std::unique_ptr<progress_meter::ProgressMeter> progress_meter;
			if (progress) {
				progress_meter = std::make_unique<progress_meter::ProgressMeter>(
						total_term_updates, "[odgi::path_linear_sgd] 1D path-guided SGD:");
			}
			using namespace std::chrono_literals; // for timing stuff
			atomic<bool> snapshot_in_progress;
			snapshot_in_progress.store(false);
			std::vector<atomic<bool>> snapshot_progress(iter_max);
			// we will produce one less snapshot compared to iterations
			snapshot_progress[0].store(true);
			bool at_least_one_path_with_more_than_one_step = false;

			for (auto &path : path_sgd_use_paths) {
#ifdef debug_path_sgd
				std::string path_name = graph.get_path_name(path);
                std::cerr << path_name << std::endl;
                std::cerr << as_integer(path) << std::endl;
#endif
#ifdef debug_path_sgd
				size_t path_len = path_index.get_path_length(path);
                std::cerr << path_name << " has length: " << path_len << std::endl;
#endif
				if (path_index.get_path_step_count(path) > 1){
					at_least_one_path_with_more_than_one_step = true;
					break;
				}
			}

			if (at_least_one_path_with_more_than_one_step){
				double w_min = (double) 1.0 / (double) (eta_max);

#ifdef debug_path_sgd
				std::cerr << "w_min " << w_min << std::endl;
#endif
				double w_max = 1.0;
				// get our schedule
				if (progress) {
					std::cerr << "[odgi::path_linear_sgd] calculating linear SGD schedule (" << w_min << " " << w_max << " "
							  << iter_max << " " << iter_with_max_learning_rate << " " << eps << ")" << std::endl;
				}
				std::vector<double> etas = path_linear_sgd_schedule_ssi(w_min,
																	w_max,
																	iter_max,
																	iter_with_max_learning_rate,
																	eps);

				// cache zipf zetas for our full path space (heavy, but one-off)
				if (progress) {
					std::cerr << "[odgi::path_linear_sgd] calculating zetas for " << (space <= space_max ? space : space_max + (space - space_max) / space_quantization_step + 1) << " zipf distributions" << std::endl;
				}

				std::vector<double> zetas((space <= space_max ? space : space_max + (space - space_max) / space_quantization_step + 1)+1);
				uint64_t last_quantized_i = 0;
#pragma omp parallel for schedule(static,1)
				for (uint64_t i = 1; i < space+1; ++i) {
					uint64_t quantized_i = i;
					uint64_t compressed_space = i;
					if (i > space_max){
						quantized_i = space_max + (i - space_max) / space_quantization_step + 1;
						compressed_space = space_max + ((i - space_max) / space_quantization_step) * space_quantization_step;
					}

					if (quantized_i != last_quantized_i){
						dirtyzipf::dirty_zipfian_int_distribution<uint64_t>::param_type z_p(1, compressed_space, theta);
						zetas[quantized_i] = z_p.zeta();

						last_quantized_i = quantized_i;
					}
				}

				// how many term updates we make
				std::atomic<uint64_t> term_updates;
				term_updates.store(0);
				// learning rate
				std::atomic<double> eta;
				eta.store(etas.front());
				// adaptive zip theta
				std::atomic<double> adj_theta;
				adj_theta.store(theta);
				// if we're in a final cooling phase (last 10%) of iterations
				std::atomic<bool> cooling;
				cooling.store(false);
				// our max delta
				std::atomic<double> Delta_max;
				Delta_max.store(0);
				// should we keep working?
				std::atomic<bool> work_todo;
				work_todo.store(true);
				// approximately what iteration we're on
				uint64_t iteration = 0;
				// launch a thread to update the learning rate, count iterations, and decide when to stop
				auto checker_lambda =
						[&]() {
							while (work_todo.load()) {
								if (term_updates.load() > min_term_updates) {
									if (snapshot) {
										if (snapshot_progress[iteration].load() || iteration == iter_max) {
											iteration++;
											if (iteration == iter_max) {
												snapshot_in_progress.store(false);
											} else {
												snapshot_in_progress.store(true);
											}
										} else {
											snapshot_in_progress.store(true);
											continue;
										}
									} else {
										iteration++;
										snapshot_in_progress.store(false);
									}
									if (iteration > iter_max) {
										work_todo.store(false);
									} else if (Delta_max.load() <= delta) { // nb: this will also break at 0
										if (progress) {
											std::cerr << "[odgi::path_linear_sgd] delta_max: " << Delta_max.load()
													  << " <= delta: "
													  << delta << ". Threshold reached, therefore ending iterations."
													  << std::endl;
										}
										work_todo.store(false);
									} else {
										eta.store(etas[iteration]); // update our learning rate
										Delta_max.store(delta); // set our delta max to the threshold
										if (iteration > first_cooling_iteration) {
											adj_theta.store(0.001);
											cooling.store(true);
										}
									}
									term_updates.store(0);
								}
								std::this_thread::sleep_for(1ms);
							}
						};

				auto worker_lambda =
						[&](uint64_t tid) {
							// everyone tries to seed with their own random data
							const std::uint64_t seed = 9399220 + tid;
							XoshiroCpp::Xoshiro256Plus gen(seed); // a nice, fast PRNG
							// FIXME these should be replaced by the bv and ssi
							// some references to literal bitvectors in the path index hmmm
							// const sdsl::bit_vector &np_bv = path_index.get_np_bv();
							const sdsl::int_vector<> &nr_iv = path_index.get_nr_iv();
							const sdsl::int_vector<> &npi_iv = path_index.get_npi_iv();
							// we'll sample from all path steps
							// we have a larger search space because of ns_bv + node_count
							std::uniform_int_distribution<uint64_t> dis_step = std::uniform_int_distribution<uint64_t>(0, ns_bv.size() - 1);
							// std::uniform_int_distribution<uint64_t> dis_step = std::uniform_int_distribution<uint64_t>(0, np_bv.size() - 1);
							std::uniform_int_distribution<uint64_t> flip(0, 1);
							uint64_t term_updates_local = 0;
							while (work_todo.load()) {
								if (!snapshot_in_progress.load()) {
									// sample the first node from all the nodes in the graph
									// pick a random position from all paths
									uint64_t node_step_index = dis_step(gen);
#ifdef debug_sample_from_nodes
									std::cerr << "step_index: " << node_step_index << std::endl;
#endif
									// if ns_bv[step_index] == 1 we have to continue, because we hit a node!
									if (ns_bv[node_step_index] == 1) {
										continue;
									}
									uint64_t handle_rank_of_step_index = ns_bv_rank(node_step_index); // because we have a compacted graph, this is the actual node id!
									uint64_t step_index = node_step_index - handle_rank_of_step_index;
									uint64_t step_rank_on_node_given_step_index = node_step_index - ns_bv_select(handle_rank_of_step_index) - 1;
									// I want to have the handle
									handle_t term_i_ssi = graph.get_handle(handle_rank_of_step_index, handle_orientation[step_index]);
									// I want to have the step
									uint64_t step_rank_on_handle = 0;
									step_handle_t step_a_ssi, step_b_ssi;
									graph.for_each_step_on_handle(
											term_i_ssi,
											[&](const step_handle_t& step) {
												if (step_rank_on_handle == step_rank_on_node_given_step_index) {
													step_a_ssi = step;
													// FIXME once we found this, step out of it
													// I am not sure what is happening here, but this return seems to abort the mission too early
													// return;
												}
												step_rank_on_handle++;
											});

									uint64_t path_i = npi_iv[step_index];
									// FIXME delete this later
									// path_handle_t path = as_path_handle(path_i);
									path_handle_t path = graph.get_path_handle_of_step(step_a_ssi);

									size_t path_step_count = graph.get_step_count(path);
									if (path_step_count == 1){
										continue;
									}
#ifdef debug_sample_from_nodes
									std::cerr << "path integer: " << path_i << std::endl;
#endif
									step_handle_t step_a, step_b;
									as_integers(step_a)[0] = path_i; // path index
									size_t s_rank = nr_iv[step_index] - 1; // step rank in path
									as_integers(step_a)[1] = s_rank;
#ifdef debug_sample_from_nodes
									std::cerr << "step rank in path: " << nr_iv[step_index]  << std::endl;
#endif
									bool cool = false;
									if (cooling.load() || flip(gen)) {
										cool = true;
										auto _theta = adj_theta.load();

										/// IN PROGRESS
										bool path_end = false;
										const uint64_t path_steps = sampled_step_index.get_path_len(path);
										uint64_t local_space;
										// FIXME we need to limit the jump length by the given path length?
										// FIXME I think the problem is how the zetas are created?!
										if (path_steps > space_max){
											local_space = space_max + (path_steps - space_max) / space_quantization_step + 1;
										}
										//dirtyzipf::dirty_zipfian_int_distribution<uint64_t>::param_type z_p(1, path_steps, _theta, zetas[local_space]);
										dirtyzipf::dirty_zipfian_int_distribution<uint64_t>::param_type z_p(1, path_steps/1000, _theta);
										dirtyzipf::dirty_zipfian_int_distribution<uint64_t> z(z_p);
										uint64_t z_i = z(gen);

										uint64_t steps_travelled = 0;
										step_handle_t last_step_in_path = graph.path_back(path);
										step_handle_t first_step_in_path = graph.path_begin(path);
										step_handle_t cur_step = step_a_ssi;

										if (flip(gen)) {
											// go backwards
											while (!path_end) {
												// did we reach the left end already?
												if ((as_integers(cur_step)[0] == as_integers(first_step_in_path)[0]) &&
													(as_integers(cur_step)[1] == as_integers(first_step_in_path)[1])) {
													cur_step = graph.path_back(path);
												} else {
													cur_step = graph.get_previous_step(cur_step);
												}
												steps_travelled++;
												// this assumes z_i is at least 1;
												if (z_i == steps_travelled) {
													path_end = true;
												}
											}
											// go forward
										} else {
											while (!path_end) {
												// did we reach the right already?
												if ((as_integers(cur_step)[0] == as_integers(last_step_in_path)[0]) &&
													(as_integers(cur_step)[1] == as_integers(last_step_in_path)[1])) {
													cur_step = graph.path_begin(path);
												} else {
													cur_step = graph.get_next_step(cur_step);
												}
												steps_travelled++;
												if (z_i == steps_travelled) {
													path_end = true;
												}
											}
										}
										step_b_ssi = cur_step;
										/// IN PROGRESS END

										// if I understand correctly, we can replace [s_rank == path_step_count - 1] with graph.has_next_step(step_a)?
										// BUT TRY TO COMPARE AGAINST THE LAST STEP, BECAUSE HAS NEXT_STEP IS SLOW
										if (s_rank > 0 && flip(gen) || s_rank == path_step_count-1) {
											// FIXME why do we check if the step_rank is larger 0?
											// go backward
											uint64_t jump_space = std::min(space, s_rank);
											uint64_t space = jump_space;
											if (jump_space > space_max){
												space = space_max + (jump_space - space_max) / space_quantization_step + 1;
											}
											dirtyzipf::dirty_zipfian_int_distribution<uint64_t>::param_type z_p(1, jump_space, _theta, zetas[space]);
											dirtyzipf::dirty_zipfian_int_distribution<uint64_t> z(z_p);
											uint64_t z_i = z(gen);
											//assert(z_i <= path_space);
											// FIXME travel there :)
											as_integers(step_b)[0] = as_integer(path);
											as_integers(step_b)[1] = s_rank - z_i;
										} else {
											// go forward
											uint64_t jump_space = std::min(space, path_step_count - s_rank - 1);
											uint64_t space = jump_space;
											if (jump_space > space_max){
												space = space_max + (jump_space - space_max) / space_quantization_step + 1;
											}
											dirtyzipf::dirty_zipfian_int_distribution<uint64_t>::param_type z_p(1, jump_space, _theta, zetas[space]);
											dirtyzipf::dirty_zipfian_int_distribution<uint64_t> z(z_p);
											uint64_t z_i = z(gen);
											//assert(z_i <= path_space);
											as_integers(step_b)[0] = as_integer(path);
											as_integers(step_b)[1] = s_rank + z_i;
										}
									} else {
										// sample randomly across the path
										std::uniform_int_distribution<uint64_t> rando(0, graph.get_step_count(path)-1);
										// FIXME very inefficient solution! Is there a more random way?!
										uint64_t step_rank_b = rando(gen);
										uint64_t n = 0;
										step_handle_t cur_step = graph.path_begin(path);
										bool lost = true;
										// has_next_step is very, very costly!
										while (lost) {
											if (n == step_rank_b) {
												cur_step = graph.get_next_step(cur_step);
												lost = false;
											}
											n++;
										}
										step_b_ssi = cur_step;

										as_integers(step_b)[0] = as_integer(path);
										as_integers(step_b)[1] = step_rank_b;
									}

									// and the graph handles, which we need to record the update
									// FIXME

									/// xp:
									//as_integers(step)[0] // path_id
									//as_integers(step)[1] // step_rank?
									/// ODGI:
									//as_integers(step)[0] // handle_rank
									//as_integers(step)[1] // path_id
									// handle_t term_i = path_index.get_handle_of_step(step_a);

									handle_t term_i = term_i_ssi;
									handle_t term_j;
									if (cool) {
										term_j = path_index.get_handle_of_step(step_b);
										// FIXME laters!
									} else {
										term_j = graph.get_handle_of_step(step_b_ssi);
									}
									term_j = graph.get_handle_of_step(step_b_ssi);

									bool update_term_i = true;
									bool update_term_j = true;


									// Check which terms we actually have to update
									if (target_sorting) {
										if (target_nodes[graph.get_id(term_i) - 1]) {
											update_term_i = false;
										}
										if (target_nodes[graph.get_id(term_j) - 1]) {
											update_term_j = false;
										}
									}
									if (!update_term_j && !update_term_i) {
										// we also have to update the number of terms here, because else we will over sample and the sorting will take much longer
										term_updates_local++;
										if (progress) {
											progress_meter->increment(1);
										}
										continue;
									}

									// adjust the positions to the node starts
									/// FIXME remove this later
									// size_t pos_in_path_a = path_index.get_position_of_step(step_a);
									size_t pos_in_path_a = sampled_step_index.get_position(step_a_ssi, graph);
									/// FIXME we use ODGI to run along until we reach the step with the rank of b
									size_t pos_in_path_b = sampled_step_index.get_position(step_b_ssi, graph);
									/*
									if (cool) {
										pos_in_path_b = path_index.get_position_of_step(step_b);
									} else {
										pos_in_path_b = sampled_step_index.get_position(step_b_ssi, graph);
									}
									*/
#ifdef debug_path_sgd
									std::cerr << "1. pos in path " << pos_in_path_a << " " << pos_in_path_b << std::endl;
#endif
									// assert(pos_in_path_a < path_index.get_path_length(path));
									// assert(pos_in_path_b < path_index.get_path_length(path));

#ifdef debug_path_sgd
									std::cerr << "2. pos in path " << pos_in_path_a << " " << pos_in_path_b << std::endl;
#endif
									// establish the term distance
									double term_dist = std::abs(
											static_cast<double>(pos_in_path_a) - static_cast<double>(pos_in_path_b));

									if (term_dist == 0) {
										continue;
										// term_dist = 1e-9;
									}
#ifdef eval_path_sgd
									std::string path_name = graph.get_path_name(path);
                                	std::cerr << path_name << "\t" << pos_in_path_a << "\t" << pos_in_path_b << "\t" << term_dist << std::endl;
#endif
									// assert(term_dist == zipf_int);
#ifdef debug_path_sgd
									std::cerr << "term_dist: " << term_dist << std::endl;
#endif
									double term_weight = 1.0 / term_dist;

									double w_ij = term_weight;
#ifdef debug_path_sgd
									std::cerr << "w_ij = " << w_ij << std::endl;
#endif
									double mu = eta.load() * w_ij;
									if (mu > 1) {
										mu = 1;
									}
									// actual distance in graph
									double d_ij = term_dist;
									// identities
									uint64_t i = number_bool_packing::unpack_number(term_i);
									uint64_t j = number_bool_packing::unpack_number(term_j);
#ifdef debug_path_sgd
									#pragma omp critical (cerr)
                                std::cerr << "nodes are " << graph.get_id(term_i) << " and " << graph.get_id(term_j) << std::endl;
#endif
									// distance == magnitude in our 1D situation
									double dx = X[i].load() - X[j].load();
									if (dx == 0) {
										dx = 1e-9; // avoid nan
									}
#ifdef debug_path_sgd
									#pragma omp critical (cerr)
                                std::cerr << "distance is " << dx << " but should be " << d_ij << std::endl;
#endif
									//double mag = dx; //sqrt(dx*dx + dy*dy);
									double mag = std::abs(dx);
#ifdef debug_path_sgd
									std::cerr << "mu " << mu << " mag " << mag << " d_ij " << d_ij << std::endl;
#endif
									// check distances for early stopping
									double Delta = mu * (mag - d_ij) / 2;
									// try until we succeed. risky.
									double Delta_abs = std::abs(Delta);
#ifdef debug_path_sgd
									#pragma omp critical (cerr)
                                std::cerr << "Delta_abs " << Delta_abs << std::endl;
#endif
									while (Delta_abs > Delta_max.load()) {
										Delta_max.store(Delta_abs);
									}
									// calculate update
									double r = Delta / mag;
									double r_x = r * dx;
#ifdef debug_path_sgd
									#pragma omp critical (cerr)
                                std::cerr << "r_x is " << r_x << std::endl;
#endif
									// update our positions (atomically)
#ifdef debug_path_sgd
									std::cerr << "before X[i] " << X[i].load() << " X[j] " << X[j].load() << std::endl;
#endif
									if (update_term_i) {
										X[i].store(X[i].load() - r_x);
									}
									if (update_term_j) {
										X[j].store(X[j].load() + r_x);
									}
#ifdef debug_path_sgd
									std::cerr << "after X[i] " << X[i].load() << " X[j] " << X[j].load() << std::endl;
#endif
									term_updates_local++;
									if (term_updates_local >= 1000) {
										term_updates += term_updates_local;
										term_updates_local = 0;
									}
									if (progress) {
										progress_meter->increment(1);
									}
								}
							}
						};

				auto snapshot_lambda =
						[&](void) {
							uint64_t iter = 0;
							while (snapshot && work_todo.load()) {
								if ((iter < iteration) && iteration != iter_max) {
									//snapshot_in_progress.store(true); // will be released again by the snapshot thread
									std::cerr << "[odgi::path_linear_sgd] snapshot thread: Taking snapshot!" << std::endl;
									// create temp file
									std::string snapshot_tmp_file = xp::temp_file::create("snapshot");
									// write to temp file
									ofstream snapshot_stream;
									snapshot_stream.open(snapshot_tmp_file);
									for (auto &x : X) {
										snapshot_stream << x << std::endl;
									}
									// push back the name of the temp file
									snapshots.push_back(snapshot_tmp_file);
									iter = iteration;
									// std::cerr << "ITER: " << iter << std::endl;
									snapshot_in_progress.store(false);
									snapshot_progress[iter].store(true);
								}
								std::this_thread::sleep_for(1ms);
							}

						};

				std::thread checker(checker_lambda);
				std::thread snapshot_thread(snapshot_lambda);

				std::vector<std::thread> workers;
				workers.reserve(nthreads);
				for (uint64_t t = 0; t < nthreads; ++t) {
					workers.emplace_back(worker_lambda, t);
				}

				for (uint64_t t = 0; t < nthreads; ++t) {
					workers[t].join();
				}

				snapshot_thread.join();

				checker.join();
			}

			if (progress) {
				progress_meter->finish();
			}

			// drop out of atomic stuff... maybe not the best way to do this
			std::vector<double> X_final(X.size());
			uint64_t i = 0;
			for (auto &x : X) {
				X_final[i++] = x.load();
			}
			return X_final;
		}

		std::vector<double> path_linear_sgd_schedule_ssi(const double &w_min,
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

		std::vector<handle_t> path_linear_sgd_order_ssi(graph_t &graph,
													const algorithms::step_index_t &sampled_step_index,
													const std::vector<path_handle_t> &path_sgd_use_paths,
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
													const uint64_t &sum_path_step_count) {
			std::vector<string> snapshots;
			std::vector<double> layout = path_linear_sgd_ssi(graph,
														 sampled_step_index,
														 path_sgd_use_paths,
														 iter_max,
														 iter_with_max_learning_rate,
														 min_term_updates,
														 delta,
														 eps,
														 eta_max,
														 theta,
														 space,
														 space_max,
														 space_quantization_step,
														 cooling_start,
														 nthreads,
														 progress,
														 snapshot,
														 snapshots,
														 target_sorting,
														 target_nodes,
														 sum_path_step_count);
			// TODO move the following into its own function that we can reuse
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
			std::vector<uint64_t> weak_components_map;
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
			// generate and write snapshot graphs
			if (snapshot) {
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
			// from layout to order
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
			std::vector<handle_t> order;
			order.reserve(graph.get_node_count());
			for (auto &layout_handle : handle_layout) {
				order.push_back(layout_handle.handle);
			}
			return order;
		}
	}
}
