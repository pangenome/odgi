#include "linear_sgd.hpp"

namespace odgi {
namespace algorithms {

std::vector<double> linear_sgd(const HandleGraph& graph,
                               const uint64_t& bandwidth,
                               const uint64_t& t_max,
                               const double& eps,
                               const double& delta,
                               const uint64_t& nthreads) {
    using namespace std::chrono_literals; // for timing stuff
     // our positions in 1D
    std::vector<std::atomic<double>> X(graph.get_node_count());
    // seed them with the graph order
    uint64_t len = 0;
    graph.for_each_handle(
        [&X,&graph,&len](const handle_t& handle) {
            // nb: we assume that the graph provides a compact handle set
            X[number_bool_packing::unpack_number(handle)].store(len);
            len += graph.get_length(handle);
        });
    // run banded BFSs in the graph to build our terms
    std::vector<sgd_term_t> terms = linear_sgd_search(graph, bandwidth);
    // locks for each term
    std::vector<std::mutex> term_mutexes(terms.size());
    // get our schedule
    std::vector<double> etas = linear_sgd_schedule(terms, t_max, eps);
    // iterate through step sizes
    /*
    for (auto &t : terms) {
        std::cerr << "considering " << graph.get_id(t.i) << " and " << graph.get_id(t.j) << " w = " << t.w << " d = " << t.d << std::endl;
    }
    */
    // how many term updates we make
    std::atomic<uint64_t> term_updates; term_updates.store(0);
    // learning rate
    std::atomic<double> eta; eta.store(etas.front());
    // our max delta
    std::atomic<double> Delta_max; Delta_max.store(0);
    // should we keep working?
    std::atomic<bool> work_todo; work_todo.store(true);
    // approximately what iteration we're on
    uint64_t iteration = 0;
    // launch a thread to update the learning rate, count iterations, and decide when to stop
    auto checker_lambda =
        [&](void) {
            while (work_todo.load()) {
                if (term_updates.load() > terms.size()) {
                    std::cerr << iteration
                              << ", eta: " << eta.load()
                              << ", Delta: " << Delta_max.load()
                              << " terms " << terms.size()
                              << " updates " << term_updates.load() << std::endl;
                    if (++iteration > t_max) {
                        work_todo.store(false);
                    } else if (Delta_max.load() <= delta) { // nb: this will also break at 0
                        work_todo.store(false);
                    } else {
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
            std::uniform_int_distribution<uint64_t> dis(0, terms.size()-1);
            while (work_todo.load()) {
                // pick a random term
                uint64_t idx = dis(gen);
                // try to acquire a mutex lock on it, and try again if we fail
                if (term_mutexes[idx].try_lock()) {
                    auto& t = terms[idx];
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
                    // distance == magnitude in our 1D situation
                    double dx = X[i].load()-X[j].load(); //, dy = X[i*2+1]-X[j*2+1];
                    //double mag = dx; //sqrt(dx*dx + dy*dy);
                    double mag = sqrt(dx*dx);
                    // check distances for early stopping
                    double Delta = mu * (mag-d_ij) / 2;
                    // try until we succeed. risky.
                    while (Delta > Delta_max.load()) {
                        Delta_max.store(Delta);
                    }
                    // calculate update
                    double r = Delta / mag;
                    double r_x = r * dx;
                    // update our positions (atomically)
                    X[i].store(X[i].load()-r_x);
                    X[j].store(X[j].load()+r_x); 
                    term_updates.store(term_updates.load()+1);
                    term_mutexes[idx].unlock();
                }
            }
        };

    std::thread checker(checker_lambda);

    std::vector<std::thread> workers; workers.reserve(nthreads);
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
    for (auto& x : X) {
        X_final[i++] = x.load();
    }
    return X_final;
}

// find pairs of handles to operate on, searching up to bandwidth steps, recording their graph distance
std::vector<sgd_term_t> linear_sgd_search(const HandleGraph& graph,
                                          const uint64_t& bandwidth) {
    std::vector<sgd_term_t> terms;
    ska::flat_hash_set<std::pair<handle_t, handle_t>> seen_pairs;
    graph.for_each_handle(
        [&](const handle_t& h) {
            uint64_t seen_steps = 0;
            ska::flat_hash_set<handle_t> seen;
            bfs(graph,
                [&](const handle_t& n, const uint64_t& depth, const uint64_t& dist) {
                    seen.insert(n);
                    if (n != h) {
                        double weight = 1.0 / ((double)dist*dist);
                        if (!seen_pairs.count(std::make_pair(h, n))
                            && !seen_pairs.count(std::make_pair(n, h))) {
                            terms.push_back(sgd_term_t(h, n, dist, weight));
                            seen_pairs.insert(std::make_pair(h, n));
                        }
                        ++seen_steps;
                    }
                },
                [&](const handle_t& n) { return seen.count(n) > 0; },
                [&](void) { return seen_steps > bandwidth; },
                { h },
                { },
                false);
        });
    // todo take unique terms
    /*
    std::sort(terms.begin(), terms.end(),
              [](const sgd_term_t& a,
                 const sgd_term_t& b) {
                  auto& a_i = as_integer(a.i);
                  auto& a_j = as_integer(a.j);
                  auto& b_i = as_integer(b.i);
                  auto& b_j = as_integer(b.j);
                  return a_i < b_i || a_i == b_i && a_j < b_j;
              });
    std::sort(terms.begin(), terms.end(),
    */
    return terms;
}


std::vector<double> linear_sgd_schedule(const std::vector<sgd_term_t> &terms,
                                        const uint64_t& t_max,
                                        const double& eps) {
    double w_min = std::numeric_limits<double>::max();
    double w_max = std::numeric_limits<double>::min();
    for (auto& term : terms) {
        auto& w = term.w;
        if (w < w_min) w_min = w;
        if (w > w_max) w_max = w;
    }
    double eta_max = 1.0 / w_min;
    double eta_min = eps / w_max;
    double lambda = log(eta_max/eta_min) / ((double)t_max-1);
    // initialize step sizes
    std::vector<double> etas;
    etas.reserve(t_max);
    for (uint64_t t=0; t<t_max; t++) {
        etas.push_back(eta_max * exp(-lambda * t));
    }
    return etas;
}

std::vector<handle_t> linear_sgd_order(const HandleGraph& graph,
                                       const uint64_t& bandwidth,
                                       const uint64_t& t_max,
                                       const double& eps,
                                       const double& delta,
                                       const uint64_t& nthreads) {
    std::vector<double> layout = linear_sgd(graph, bandwidth, t_max, eps, delta, nthreads);
    std::vector<std::pair<double, handle_t>> layout_handles;
    uint64_t i = 0;
    graph.for_each_handle([&i,&layout,&layout_handles](const handle_t& handle) {
                              layout_handles.push_back(
                                  std::make_pair(
                                      layout[i++],
                                      handle));
                          });
    std::sort(layout_handles.begin(), layout_handles.end(),
              [&](const std::pair<double, handle_t>& a,
                  const std::pair<double, handle_t>& b) {
                  return a.first < b.first
                                   || (a.first == b.first
                                       && as_integer(a.second) < as_integer(b.second));
              });
    std::vector<handle_t> order;
    for (auto& p : layout_handles) {
        order.push_back(p.second);
    }
    return order;
}

}
}
