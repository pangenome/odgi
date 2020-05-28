#include "path_sgd.hpp"

namespace odgi {
    namespace algorithms {

        std::vector<double> path_linear_sgd(const PathHandleGraph &graph,
                                            const xp::XP &path_index,
                                            const std::vector<std::string> path_sgd_use_paths,
                                            const uint64_t &iter_max,
                                            const uint64_t &term_updates,
                                            const double &delta,
                                            const double &eps,
                                            const double &alpha,
                                            const uint64_t &nthreads) {
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
            /*// run banded BFSs in the graph to build our terms
            std::vector<sgd_term_t> terms;
            if (use_paths) {
                terms = linear_sgd_path_search(graph, bandwidth, sampling_rate);
            } else {
                terms = linear_sgd_search(graph, bandwidth, sampling_rate);
            }
            // locks for each term
            std::vector<std::mutex> term_mutexes(terms.size());
            // get our schedule
            std::vector<double> etas = linear_sgd_schedule(terms, t_max, eps);
            // iterate through step sizes
            *//*
            for (auto &t : terms) {
                std::cerr << "considering " << graph.get_id(t.i) << " and " << graph.get_id(t.j) << " w = " << t.w << " d = " << t.d << std::endl;
            }
            *//*
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
                            if (term_updates.load() > terms.size()) {
                                if (++iteration > t_max) {
                                    work_todo.store(false);
                                } else if (Delta_max.load() <= delta) { // nb: this will also break at 0
                                    work_todo.store(false);
                                } else {
//#pragma omp critical (cerr)
                                    *//*
                                    std::cerr << iteration
                                              << ", eta: " << eta.load()
                                              << ", Delta: " << Delta_max.load()
                                              << " terms " << terms.size()
                                              << " updates " << term_updates.load() << std::endl;
                                    *//*
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
                        std::random_device rd;
                        std::mt19937 gen(rd());
                        std::uniform_int_distribution<uint64_t> dis(0, terms.size() - 1);
                        while (work_todo.load()) {
                            // pick a random term
                            uint64_t idx = dis(gen);
                            // try to acquire a mutex lock on it, and try again if we fail
                            if (term_mutexes[idx].try_lock()) {
                                auto &t = terms[idx];
                                double w_ij = t.w;
                                double mu = eta.load() * w_ij;
                                if (mu > 1) {
                                    mu = 1;
                                }
                                // actual distance in graph
                                double d_ij = t.d;
                                // identities
                                uint64_t i = number_bool_packing::unpack_number(t.i);
                                uint64_t j = number_bool_packing::unpack_number(t.j);
//#pragma omp critical (cerr)
//                  std::cerr << "nodes are " << graph.get_id(t.i) << " and " << graph.get_id(t.j) << std::endl;
                                // distance == magnitude in our 1D situation
                                double dx = X[i].load() - X[j].load();
                                if (dx == 0) {
                                    dx = 1e-9; // avoid nan
                                }
//#pragma omp critical (cerr)
//                    std::cerr << "distance is " << dx << " but should be " << d_ij << std::endl;
                                //double mag = dx; //sqrt(dx*dx + dy*dy);
                                double mag = sqrt(dx * dx);
                                // check distances for early stopping
                                double Delta = mu * (mag - d_ij) / 2;
                                // try until we succeed. risky.
                                double Delta_abs = std::abs(Delta);
//#pragma omp critical (cerr)
//                    std::cerr << "Delta_abs " << Delta_abs << std::endl;
                                while (Delta_abs > Delta_max.load()) {
                                    Delta_max.store(Delta_abs);
                                }
                                // calculate update
                                double r = Delta / mag;
                                double r_x = r * dx;
//#pragma omp critical (cerr)
//                    std::cerr << "r_x is " << r_x << std::endl;
                                // update our positions (atomically)
                                X[i].store(X[i].load() - r_x);
                                X[j].store(X[j].load() + r_x);
                                term_updates.store(term_updates.load() + 1);
                                term_mutexes[idx].unlock();
                            }
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
*/
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
            double eta_max = 1.0 / w_min;
            double eta_min = eps / w_max;
            double lambda = log(eta_max / eta_min) / ((double) iter_max - 1);
            // initialize step sizes
            std::vector<double> etas;
            etas.reserve(iter_max);
            for (uint64_t t = 0; t < iter_max; t++) {
                etas.push_back(eta_max * exp(-lambda * t));
            }
            return etas;
        }

        std::vector<handle_t> path_linear_sgd_order(const PathHandleGraph &graph,
                                               const xp::XP &path_index,
                                               const std::vector<std::string> path_sgd_use_paths,
                                               const uint64_t &iter_max,
                                               const uint64_t &term_updates,
                                               const double &delta,
                                               const double &eps,
                                               const double &alpha,
                                               const uint64_t &nthreads) {
            std::vector<double> layout = path_linear_sgd(graph,
                                                         path_index,
                                                         path_sgd_use_paths,
                                                         iter_max,
                                                         term_updates,
                                                         delta,
                                                         eps,
                                                         alpha,
                                                         nthreads);
            std::vector<std::pair<double, handle_t>> layout_handles;
            uint64_t i = 0;
            graph.for_each_handle([&i, &layout, &layout_handles](const handle_t &handle) {
                layout_handles.push_back(
                        std::make_pair(
                                layout[i++],
                                handle));
            });
            std::sort(layout_handles.begin(), layout_handles.end(),
                      [&](const std::pair<double, handle_t> &a,
                          const std::pair<double, handle_t> &b) {
                          return a.first < b.first
                                 || (a.first == b.first
                                     && as_integer(a.second) < as_integer(b.second));
                      });
            std::vector<handle_t> order;
            for (auto &p : layout_handles) {
                order.push_back(p.second);
            }
            return order;
        }

    }
}
