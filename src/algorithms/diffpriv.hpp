#pragma once

#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <vector>
#include <set>
#include <deque>
#include <random>
#include <iostream>
#include <map>
#include <atomic>
#include <thread>
#include <sstream>
#include "hash_map.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

void diff_priv(
    const PathHandleGraph& graph,
    MutablePathDeletableHandleGraph& priv,
    //PathHandleGraph& priv,
    const double epsilon,
    const double target_coverage,
    const double min_haplotype_freq,
    const uint64_t bp_limit,
    const uint64_t nthreads);

}
}
