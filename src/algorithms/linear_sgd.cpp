#include "linear_sgd.hpp"

namespace odgi {
namespace algorithms {

std::vector<double> linear_sgd(const HandleGraph& graph,
                               const uint64_t& bandwidth,
                               const uint64_t& t_max,
                               const double& eps,
                               const double& delta) {
    // our positions in 1D
    std::vector<double> X(graph.get_node_count());
    // seed them with the graph order
    uint64_t len = 0;
    graph.for_each_handle(
        [&X,&graph,&len](const handle_t& handle) {
            // nb: we assume that the graph provides a compact handle set
            X[number_bool_packing::unpack_number(handle)] = len;
            len += graph.get_length(handle);
        });
    // run banded BFSs in the graph to build our terms
    std::vector<sgd_term_t> terms = linear_sgd_search(graph, bandwidth);
    // get our schedule
    std::vector<double> etas = linear_sgd_schedule(terms, t_max, eps);
    // iterate through step sizes
    /*
    for (auto &t : terms) {
        std::cerr << "considering " << graph.get_id(t.i) << " and " << graph.get_id(t.j) << " w = " << t.w << " d = " << t.d << std::endl;
    }
    */
    uint64_t iteration = 0;
    for (double eta : etas) {
        // shuffle terms
        std::random_shuffle(terms.begin(), terms.end());
        // our max update
        double Delta_max = 0;
        for (auto &t : terms) {
            //std::cerr << "considering " << graph.get_id(t.i) << " and " << graph.get_id(t.j) << " w = " << t.w << " d = " << t.d << std::endl;
            // cap step size
            double w_ij = t.w;
            double mu = eta * w_ij;
            if (mu > 1) {
                mu = 1;
            }
            // actual distance in graph
            double d_ij = t.d;
            // identities
            uint64_t i = number_bool_packing::unpack_number(t.i);
            uint64_t j = number_bool_packing::unpack_number(t.j);
            // distance == magnitude in our 1D situation
            double dx = X[i]-X[j]; //, dy = X[i*2+1]-X[j*2+1];
            //double mag = dx; //sqrt(dx*dx + dy*dy);
            double mag = sqrt(dx*dx);
            // check distances for early stopping
            double Delta = mu * (mag-d_ij) / 2;
            if (Delta > Delta_max) {
                Delta_max = Delta;
            }
            // calculate update
            double r = Delta / mag;
            double r_x = r * dx;
            // update our positions
            X[i] -= r_x;
            X[j] += r_x;
        }
        std::cerr << ++iteration << ", eta: " << eta << ", Delta: " << Delta_max << std::endl;
        if (Delta_max < delta) {
            return X;
        }
    }
    return X;
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
                                       const double& delta) {
    std::vector<double> layout = linear_sgd(graph, bandwidth, t_max, eps, delta);
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
