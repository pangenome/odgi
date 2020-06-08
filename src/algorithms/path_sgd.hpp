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
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/util.hpp>
#include "hash_map.hpp"
#include "xp.hpp"
#include "sgd_term.hpp"
#include "IntervalTree.h"
#include <zipfian_int_distribution.h>
#include <iomanip>

namespace odgi {
    namespace algorithms {

        using namespace handlegraph;

        /// use SGD driven, by path guided, and partly zipfian distribution sampled pairwise distances to obtain a 1D linear layout of the graph that respects its topology
        std::vector<double> path_linear_sgd(const PathHandleGraph &graph,
                                            const xp::XP &path_index,
                                            const std::set<std::string> path_sgd_use_paths,
                                            const uint64_t &iter_max,
                                            const uint64_t &min_term_updates,
                                            const double &delta,
                                            const double &eps,
                                            const double &theta,
                                            const uint64_t &space,
                                            const uint64_t &nthreads,
                                            const bool &progress);

        /// our learning schedule
        std::vector<double> path_linear_sgd_schedule(const double &w_min,
                                                const double &w_max,
                                                const uint64_t &iter_max,
                                                const double &eps);

        std::vector<handle_t> path_linear_sgd_order(const PathHandleGraph &graph,
                                                    const xp::XP &path_index,
                                                    const std::set<std::string> path_sgd_use_paths,
                                                    const uint64_t &iter_max,
                                                    const uint64_t &min_term_updates,
                                                    const double &delta,
                                                    const double &eps,
                                                    const double &theta,
                                                    const uint64_t &space,
                                                    const uint64_t &nthreads,
                                                    const bool &progress);

    }
}
