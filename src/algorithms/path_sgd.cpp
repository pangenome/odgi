#include "path_sgd.hpp"

// #define debug_path_sgd
// #define eval_path_sgd
// #define debug_schedule
namespace odgi {
    namespace algorithms {

        std::vector<double> path_linear_sgd(const PathHandleGraph &graph,
                                            const xp::XP &path_index,
                                            const std::vector<path_handle_t>& path_sgd_use_paths,
                                            const uint64_t &iter_max,
                                            const uint64_t &min_term_updates,
                                            const double &delta,
                                            const double &eps,
                                            const double &theta,
                                            const uint64_t &space,
                                            const uint64_t &nthreads,
                                            const bool &progress) {
#ifdef debug_path_sgd
            std::cerr << "iter_max: " << iter_max << std::endl;
            std::cerr << "min_term_updates: " << min_term_updates << std::endl;
            std::cerr << "delta: " << delta << std::endl;
            std::cerr << "eps: " << eps << std::endl;
            std::cerr << "theta: " << theta << std::endl;
            std::cerr << "space: " << space << std::endl;
#endif
            using namespace std::chrono_literals; // for timing stuff
            // our positions in 1D
            std::vector<std::atomic<double>> X(graph.get_node_count());
            // seed them with the graph order
            uint64_t len = 0;
            graph.for_each_handle(
                    [&X, &graph, &len](const handle_t &handle) {
                        // nb: we assume that the graph provides a compact handle set
                        X[number_bool_packing::unpack_number(handle)].store(len);
                        len += graph.get_length(handle);
                    });
            // the longest path length measured in nucleotides
            size_t longest_path_in_nucleotides = 0;
            // the total path length in nucleotides
            size_t total_path_len_in_nucleotides = 0;
            // here we store all path nucleotides lengths so we know later from which path we sampled our random position from
            IITree<uint64_t, path_handle_t> path_nucleotide_tree;
            // iterate over all relevant path_handles:
            //  1. build the interval tree
            //  2. find out the longest path in nucleotides and store this number size_t
            //  3. add the current path length to the total length
            for (auto& path : path_sgd_use_paths) {
#ifdef debug_path_sgd
                std::string path_name = graph.get_path_name(path);
                std::cerr << path_name << std::endl;
                std::cerr << as_integer(path) << std::endl;
#endif
                size_t path_len = path_index.get_path_length(path);
#ifdef debug_path_sgd
                std::cerr << path_name << " has length: " << path_len << std::endl;
#endif
                path_nucleotide_tree.add(total_path_len_in_nucleotides,
                                         total_path_len_in_nucleotides + path_len,
                                         path);

                if (path_len > longest_path_in_nucleotides) {
                    longest_path_in_nucleotides = path_len;
                }
                total_path_len_in_nucleotides += path_len;
            }
            path_nucleotide_tree.index();

            double w_min = (double) 1.0 / (double) (longest_path_in_nucleotides * longest_path_in_nucleotides);
#ifdef debug_path_sgd
            std::cerr << "w_min " << w_min << std::endl;
#endif
            double w_max = 1.0;
            // get our schedule
            std::vector<double> etas = path_linear_sgd_schedule(w_min, w_max, iter_max, eps);
            // initialize Zipfian distrubution so we only have to calculate zeta once
            zipfian_int_distribution<uint64_t>::param_type p(1, space, theta);
            zipfian_int_distribution<uint64_t> zipfian(p);
            // how many term updates we make
            std::atomic<uint64_t> term_updates;
            term_updates.store(0);
            // learning rate
            std::atomic<double> eta;
            eta.store(etas.front());
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
                    [&](void) {
                        while (work_todo.load()) {
                            if (term_updates.load() > min_term_updates) {
                                if (++iteration > iter_max) {
                                    work_todo.store(false);
                                } else if (Delta_max.load() <= delta) { // nb: this will also break at 0
                                    if (progress) {
                                        std::cerr << "[path sgd sort]: delta_max: " << Delta_max.load() << " <= delta: "
                                                  << delta << ". Threshold reached, therefore ending iterations."
                                                  << std::endl;
                                    }
                                    work_todo.store(false);
                                } else {
                                    if (progress) {
                                        double percent_progress = ((double) iteration / (double) iter_max) * 100.0;
                                        std::cerr << std::fixed << std::setprecision(2) << "[path sgd sort]: "
                                                  << percent_progress << "% progress: "
                                                                         "iteration: " << iteration <<
                                                  ", eta: " << eta.load() <<
                                                  ", delta_max: " << Delta_max.load() <<
                                                  ", number of updates: " << term_updates.load() << std::endl;
                                    }
                                    eta.store(etas[iteration]); // update our learning rate
                                    Delta_max.store(delta); // set our delta max to the threshold
                                }
                                term_updates.store(0);
                            }
                            std::this_thread::sleep_for(1ms);
                        }
                    };

            auto worker_lambda =
                    [&](uint64_t tid) {
                        // everyone tries to seed with their own random data
                        std::array<uint64_t, 2> seed_data = {(uint64_t)std::time(0), tid};
                        std::seed_seq sseq(std::begin(seed_data), std::end(seed_data));
                        std::mt19937_64 gen(sseq);
                        std::uniform_int_distribution<uint64_t> dis(0, total_path_len_in_nucleotides-1);
                        std::uniform_int_distribution<uint64_t> flip(0, 1);
                        while (work_todo.load()) {
                            // pick a random position from all paths
                            uint64_t pos = dis(gen);
                            // use our interval tree to get the path handle and path nucleotide position of the picked position
                            //std::vector<Interval<size_t, path_handle_t> > result;
                            //result = path_nucleotide_tree.findOverlapping(pos, pos);
                            std::vector<size_t> a;
                            path_nucleotide_tree.overlap(pos, pos+1, a);
                            if (a.empty()) {
                                std::cerr << "[odgi::path_sgd] no overlapping intervals at position " << pos << std::endl;
                                exit(1);
                            }
                            auto& p = a[0];
                            path_handle_t path = path_nucleotide_tree.data(p);
                            size_t path_start_pos = path_nucleotide_tree.start(p);
                            // size_t path_end_pos = result[0].stop;
                            size_t path_len = path_index.get_path_length(path) - 1;
                            // we have a 0-based positioning in the path index
                            size_t pos_in_path_a = pos - path_start_pos;
                            uint64_t zipf_int = zipfian(gen);
#ifdef debug_path_sgd
                            std::cerr << "random pos: " << pos << std::endl;
                            std::cerr << "path_start_pos: " << path_start_pos << std::endl;
                            std::cerr << "pos_in_path_a: " << pos_in_path_a << std::endl;
                            std::cerr << "path_len: " << path_len << std::endl;
                            //std::cerr << "zipf: " << zipf_int << std::endl;
#endif
                            size_t pos_in_path_b = pos_in_path_a;
                            if (flip(gen)) {
                                if (zipf_int > pos_in_path_a) {
                                    if (pos_in_path_a == 0) {
                                        continue;
                                    } else {
                                        zipf_int %= pos_in_path_a;
                                    }
                                }
                                pos_in_path_b -= zipf_int;
                            } else {
                                if (zipf_int > path_len - pos_in_path_a) {
                                    if (path_len - pos_in_path_a == 0) {
                                        continue;
                                    } else {
                                        zipf_int %= path_len - pos_in_path_a;
                                    }
                                }
                                pos_in_path_b += zipf_int;
                            }
#ifdef debug_path_sgd
                            std::cerr << "zipf: " << zipf_int << std::endl;
                            std::cerr << "pos_in_path_a: " << pos_in_path_a << std::endl;
                            std::cerr << "pos_in_path_b: " << pos_in_path_b << std::endl;
#endif
                            // get the step handles
                            step_handle_t step_a = path_index.get_step_at_position(path, pos_in_path_a);
                            step_handle_t step_b = path_index.get_step_at_position(path, pos_in_path_b);

                            // and the graph handles, which we need to record the update
                            handle_t term_i = path_index.get_handle_of_step(step_a);
                            handle_t term_j = path_index.get_handle_of_step(step_b);

                            // adjust the positions to the node starts
                            pos_in_path_a = path_index.get_position_of_step(step_a);
                            pos_in_path_b = path_index.get_position_of_step(step_b);
#ifdef debug_path_sgd
                            std::cerr << "1. pos in path " << pos_in_path_a << " " << pos_in_path_b << std::endl;
#endif
                            // assert(pos_in_path_a < path_index.get_path_length(path));
                            // assert(pos_in_path_b < path_index.get_path_length(path));

                            // and adjust to account for our relative orientation
                            if (graph.get_is_reverse(term_i)) {
                                pos_in_path_a += graph.get_length(term_i);
                            }
                            if (graph.get_is_reverse(term_j)) {
                                pos_in_path_b += graph.get_length(term_j);
                            }
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
                            std::string path_name = path_index.get_path_name(path);
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
                            std::cerr <
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
                            X[i].store(X[i].load() - r_x);
                            X[j].store(X[j].load() + r_x);
#ifdef debug_path_sgd
                            std::cerr << "after X[i] " << X[i].load() << " X[j] " << X[j].load() << std::endl;
#endif
                            term_updates.store(term_updates.load() + 1);
                        }
                    };

            std::thread checker(checker_lambda);

            std::vector<std::thread> workers;
            workers.reserve(nthreads);
            for (uint64_t t = 0; t < nthreads; ++t) {
                workers.emplace_back(worker_lambda, t);
            }

            for (uint64_t t = 0; t < nthreads; ++t) {
                workers[t].join();
            }
            checker.join();
            // drop out of atomic stuff... maybe not the best way to do this
            std::vector<double> X_final(X.size());
            uint64_t i = 0;
            for (auto &x : X) {
                X_final[i++] = x.load();
            }
            return X_final;
        }

        std::vector<double> path_linear_sgd_schedule(const double &w_min,
                                                     const double &w_max,
                                                     const uint64_t &iter_max,
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
            etas.reserve(iter_max);
#ifdef debug_schedule
            std::cerr << "etas: ";
#endif
            for (uint64_t t = 0; t < iter_max; t++) {
                etas.push_back(eta_max * exp(-lambda * t));
#ifdef debug_schedule
                std::cerr << etas.back() << ", ";
#endif
            }
#ifdef debug_schedule
            std::cerr << std::endl;
#endif
            return etas;
        }

        std::vector<double> deterministic_path_linear_sgd(const PathHandleGraph &graph,
                                                          const xp::XP &path_index,
                                                          const std::vector<path_handle_t>& path_sgd_use_paths,
                                                          const uint64_t &iter_max,
                                                          const uint64_t &min_term_updates,
                                                          const double &delta,
                                                          const double &eps,
                                                          const double &theta,
                                                          const uint64_t &space,
                                                          const std::string &seeding_string,
                                                          const bool &progress) {
            using namespace std::chrono_literals; // for timing stuff
            // our positions in 1D
            std::vector<std::atomic<double>> X(graph.get_node_count());
            // seed them with the graph order
            uint64_t len = 0;
            graph.for_each_handle(
                    [&X, &graph, &len](const handle_t &handle) {
                        // nb: we assume that the graph provides a compact handle set
                        X[number_bool_packing::unpack_number(handle)].store(len);
                        len += graph.get_length(handle);
                    });
            // the longest path length measured in nucleotides
            size_t longest_path_in_nucleotides = 0;
            // the total path length in nucleotides
            size_t total_path_len_in_nucleotides = 0;
            // here we store all path nucleotides lengths so we know later from which path we sampled our random position from
            IITree<uint64_t, path_handle_t> path_nucleotide_tree;
            // iterate over all relevant path_handles:
            //  1. build the interval tree
            //  2. find out the longest path in nucleotides and store this number size_t
            //  3. add the current path length to the total length
            for (auto& path : path_sgd_use_paths) {
#ifdef debug_path_sgd
                std::string path_name = graph.get_path_name(path);
                std::cerr << path_name << std::endl;
                std::cerr << as_integer(path) << std::endl;
#endif
                size_t path_len = path_index.get_path_length(path);
#ifdef debug_path_sgd
                std::cerr << path_name << " has length: " << path_len << std::endl;
#endif
                path_nucleotide_tree.add(total_path_len_in_nucleotides,
                                         total_path_len_in_nucleotides + path_len,
                                         path);

                if (path_len > longest_path_in_nucleotides) {
                    longest_path_in_nucleotides = path_len;
                }
                total_path_len_in_nucleotides += path_len;
            }
            path_nucleotide_tree.index();

            double w_min = (double) 1.0 / (double) (longest_path_in_nucleotides * longest_path_in_nucleotides);
#ifdef debug_path_sgd
            std::cerr << "w_min " << w_min << std::endl;
#endif
            double w_max = 1.0;
            // get our schedule
            std::vector<double> etas = path_linear_sgd_schedule(w_min, w_max, iter_max, eps);
            // initialize Zipfian distrubution so we only have to calculate zeta once
            zipfian_int_distribution<uint64_t>::param_type p(1, space, theta);
            zipfian_int_distribution<uint64_t> zipfian(p);
            // how many term updates we make
            std::atomic<uint64_t> term_updates;
            term_updates.store(0);
            // learning rate
            std::atomic<double> eta;
            eta.store(etas.front());
            // our max delta
            std::atomic<double> Delta_max;
            Delta_max.store(0);
            // approximately what iteration we're on
            uint64_t iteration = 0;
            // seed with the given string
            std::seed_seq seed(seeding_string.begin(), seeding_string.end());
            std::mt19937 gen(seed);
            std::uniform_int_distribution<uint64_t> dis(0, total_path_len_in_nucleotides-1);
            std::uniform_int_distribution<uint64_t> flip(0, 1);
            for (uint64_t iteration = 0; iteration < iter_max; iteration++) {
                for (uint64_t term_update = 0; term_update < min_term_updates; term_update++) {
                    // pick a random position from all paths
                    uint64_t pos = dis(gen);
#ifdef debug_path_sgd
                    std::cerr << "uniform_position: " << pos << std::endl;
#endif
                    // use our interval tree to get the path handle and path nucleotide position of the picked position
                    std::vector<size_t> a;
                    path_nucleotide_tree.overlap(pos, pos+1, a);
                    if (a.empty()) {
                        std::cerr << "[odgi::path_sgd] no overlapping intervals at position " << pos << std::endl;
                        exit(1);
                    }
                    auto& p = a[0];
                    path_handle_t path = path_nucleotide_tree.data(p);
                    size_t path_start_pos = path_nucleotide_tree.start(p);
                    // size_t path_end_pos = result[0].stop;
                    size_t path_len = path_index.get_path_length(path) - 1;
                    // we have a 0-based positioning in the path index
                    size_t pos_in_path_a = pos - path_start_pos;
                    uint64_t zipf_int = zipfian(gen);
#ifdef debug_path_sgd
                    std::cerr << "random pos: " << pos << std::endl;
                    std::cerr << "path_start_pos: " << path_start_pos << std::endl;
                    std::cerr << "pos_in_path_a: " << pos_in_path_a << std::endl;
                    std::cerr << "path_len: " << path_len << std::endl;
                    //std::cerr << "zipf: " << zipf_int << std::endl;
#endif
                    size_t pos_in_path_b = pos_in_path_a;
                    if (flip(gen)) {
                        if (zipf_int > pos_in_path_a) {
                            if (pos_in_path_a == 0) {
                                continue;
                            } else {
                                zipf_int %= pos_in_path_a;
                            }
                        }
                        pos_in_path_b -= zipf_int;
                    } else {
                        if (zipf_int > path_len - pos_in_path_a) {
                            if (path_len - pos_in_path_a == 0) {
                                continue;
                            } else {
                                zipf_int %= path_len - pos_in_path_a;
                            }
                        }
                        pos_in_path_b += zipf_int;
                    }
#ifdef debug_path_sgd
                    std::cerr << "zipf: " << zipf_int << std::endl;
                    std::cerr << "pos_in_path_a: " << pos_in_path_a << std::endl;
                    std::cerr << "pos_in_path_b: " << pos_in_path_b << std::endl;
#endif
                    // get the step handles
                    step_handle_t step_a = path_index.get_step_at_position(path, pos_in_path_a);
                    step_handle_t step_b = path_index.get_step_at_position(path, pos_in_path_b);

                    // and the graph handles, which we need to record the update
                    handle_t term_i = path_index.get_handle_of_step(step_a);
                    handle_t term_j = path_index.get_handle_of_step(step_b);

                    // adjust the positions to the node starts
                    pos_in_path_a = path_index.get_position_of_step(step_a);
                    pos_in_path_b = path_index.get_position_of_step(step_b);
#ifdef debug_path_sgd
                    std::cerr << "1. pos in path " << pos_in_path_a << " " << pos_in_path_b << std::endl;
#endif
                    // assert(pos_in_path_a < path_index.get_path_length(path));
                    // assert(pos_in_path_b < path_index.get_path_length(path));

                    // and adjust to account for our relative orientation
                    if (graph.get_is_reverse(term_i)) {
                        pos_in_path_a += graph.get_length(term_i);
                    }
                    if (graph.get_is_reverse(term_j)) {
                        pos_in_path_b += graph.get_length(term_j);
                    }
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
                    std::string path_name = path_index.get_path_name(path);
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
                    double mag = sqrt(dx * dx);
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
                    X[i].store(X[i].load() - r_x);
                    X[j].store(X[j].load() + r_x);
#ifdef debug_path_sgd
                    std::cerr << "after X[i] " << X[i].load() << " X[j] " << X[j].load() << std::endl;
#endif
                    term_updates.store(term_updates.load() + 1);
                }
                if (Delta_max.load() <= delta) {
                    if (progress) {
                        std::cerr << "[path sgd sort]: delta_max: " << Delta_max.load() << " <= delta: " << delta << ". Threshold reached, therefore ending iterations." << std::endl;
                    }
                    break;
                } else {
                    if (progress) {
                        double percent_progress = ((double) (iteration + 1) / (double) iter_max) * 100.0;
                        std::cerr << std::fixed << std::setprecision(2) << "[path sgd sort]: " << percent_progress
                                  << "% progress: "
                                     "iteration: " << (iteration + 1) <<
                                  ", eta: " << eta.load() <<
                                  ", delta: " << Delta_max.load() <<
                                  ", number of updates: " << term_updates.load() << std::endl;
                    }
                    eta.store(etas[iteration]); // update our learning rate
                    Delta_max.store(delta); // set our delta max to the threshold
                }
                term_updates.store(0);
            }

            // drop out of atomic stuff... maybe not the best way to do this
            std::vector<double> X_final(X.size());
            uint64_t i = 0;
            for (auto &x : X) {
                X_final[i++] = x.load();
            }
            return X_final;
        }

        std::vector<handle_t> path_linear_sgd_order(const PathHandleGraph &graph,
                                                    const xp::XP &path_index,
                                                    const std::vector<path_handle_t>& path_sgd_use_paths,
                                                    const uint64_t &iter_max,
                                                    const uint64_t &min_term_updates,
                                                    const double &delta,
                                                    const double &eps,
                                                    const double &theta,
                                                    const uint64_t &space,
                                                    const uint64_t &nthreads,
                                                    const bool &progress,
                                                    const std::string &seed) {
            std::vector<double> layout;
            if (nthreads == 1) {
                layout = deterministic_path_linear_sgd(graph,
                                                       path_index,
                                                       path_sgd_use_paths,
                                                       iter_max,
                                                       min_term_updates,
                                                       delta,
                                                       eps,
                                                       theta,
                                                       space,
                                                       seed,
                                                       progress);
            } else {
                layout = path_linear_sgd(graph,
                                         path_index,
                                         path_sgd_use_paths,
                                         iter_max,
                                         min_term_updates,
                                         delta,
                                         eps,
                                         theta,
                                         space,
                                         nthreads,
                                         progress);
            }
#ifdef debug_components
            std::cerr << "node count: " << graph.get_node_count() << std::endl;
#endif
            // refine order by weakly connected components
            std::vector<ska::flat_hash_set<handlegraph::nid_t>> weak_components = algorithms::weakly_connected_components(&graph);
#ifdef debug_components
            std::cerr << "components count: " << weak_components.size() << std::endl;
#endif
            std::vector<std::pair<double, uint64_t>> weak_component_order;
            for (int i = 0; i < weak_components.size(); i++) {
                auto& weak_component = weak_components[i];
                uint64_t id_sum = 0;
                for (auto node_id : weak_component) {
                    id_sum += node_id;
                }
                double avg_id = id_sum / (double)weak_component.size();
                weak_component_order.push_back(std::make_pair(avg_id, i));
            }
            std::sort(weak_component_order.begin(), weak_component_order.end());
            std::vector<uint64_t> weak_component_id; // maps rank to "id" based on the orignial sorted order
            weak_component_id.resize(weak_component_order.size());
            uint64_t component_id = 0;
            for (auto& component_order : weak_component_order) {
                weak_component_id[component_order.second] = component_id++;
            }
            std::vector<uint64_t> weak_components_map;
            weak_components_map.resize(graph.get_node_count());
            // reserve the space we need
            for (int i = 0; i < weak_components.size(); i++) {
                auto& weak_component = weak_components[i];
                // store for each node identifier to component start index
                for (auto node_id : weak_component) {
                    weak_components_map[node_id-1] = weak_component_id[i];
                }
#ifdef debug_components
                std::cerr << "weak_component.size(): " << weak_component.size() << std::endl;
                std::cerr << "component_index: " << i << std::endl;
#endif
            }
            weak_components_map.clear();
            std::vector<handle_layout_t> handle_layout;
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
                      [&](const handle_layout_t& a,
                          const handle_layout_t& b) {
                          return a.weak_component < b.weak_component
                                   || (a.weak_component == b.weak_component
                                       && a.pos < b.pos
                                       || (a.pos == b.pos
                                           && as_integer(a.handle) < as_integer(b.handle)));
                      });
            std::vector<handle_t> order;
            order.reserve(graph.get_node_count());
            for (auto& layout_handle : handle_layout) {
                order.push_back(layout_handle.handle);
            }
            return order;
        }

    }
}
