#pragma once

#include <vector>
#include <string>
#include <cassert>
#include <algorithm>
#include <random>
#include <set>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/util.hpp>
#include "hash_map.hpp"
#include "bfs.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

struct sgd_term_t {
    handle_t i, j;
    double d, w;
    sgd_term_t(const handle_t& i, const handle_t& j, const double& d, const double& w) : i(i), j(j), d(d), w(w) {}
};

// use SGD driven by banded pairwise distances to obtain a linear layout of the graph that respects its topology
std::vector<double> linear_sgd(const HandleGraph& graph,
                               const uint64_t& bandwidth,
                               const uint64_t& t_max,
                               const double& eps,
                               const double& delta,
                               const uint64_t& nthreads);

// find pairs of handles to operate on, searching up to bandwidth steps, recording their graph distance
std::vector<sgd_term_t> linear_sgd_search(const HandleGraph& graph, const uint64_t& bandwidth);

// our learning schedule
std::vector<double> linear_sgd_schedule(const std::vector<sgd_term_t>& terms,
                                        const uint64_t& t_max,
                                        const double& eps);

std::vector<handle_t> linear_sgd_order(const HandleGraph& graph,
                                       const uint64_t& bandwidth,
                                       const uint64_t& t_max,
                                       const double& eps,
                                       const double& delta,
                                       const uint64_t& nthreads);

}
}
