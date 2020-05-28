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
#include "sgd_term.hpp"
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <bf/all.hpp>
#include "hash_map.hpp"
#include "bfs.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

// use SGD driven by banded pairwise distances to obtain a linear layout of the graph that respects its topology
std::vector<double> linear_sgd(const PathHandleGraph& graph,
                               const uint64_t& bandwidth,
                               const double& sampling_rate,
                               const bool& use_paths,
                               const uint64_t& t_max,
                               const double& eps,
                               const double& delta,
                               const uint64_t& nthreads);


// find pairs of handles to operate on, searching up to bandwidth steps, recording their graph distance
std::vector<sgd_term_t> linear_sgd_search(const HandleGraph& graph,
                                          const uint64_t& bandwidth,
                                          const double& sampling_rate);

// find pairs of handles to operate on using path iteration to establish distances
std::vector<sgd_term_t> linear_sgd_path_search(const PathHandleGraph& graph,
                                               const uint64_t& bandwidth,
                                               const double& sampling_rate);

// our learning schedule
std::vector<double> linear_sgd_schedule(const std::vector<sgd_term_t>& terms,
                                        const uint64_t& t_max,
                                        const double& eps);

std::vector<handle_t> linear_sgd_order(const PathHandleGraph& graph,
                                       const uint64_t& bandwidth,
                                       const double& sampling_rate,
                                       const bool& use_paths,
                                       const uint64_t& t_max,
                                       const double& eps,
                                       const double& delta,
                                       const uint64_t& nthreads);

}
}
