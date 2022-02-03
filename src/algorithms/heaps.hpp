#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <random>
#include <omp.h>
#include "hash_map.hpp"
#include "position.hpp"
#include <handlegraph/types.hpp>
#include <handlegraph/iteratee.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>

namespace odgi {

using namespace handlegraph;

namespace algorithms {

/// For each permutation of the path groups
/// we call func with a vector that is the fraction of the pangenome covered when we've considered N groups in the permutation
void for_each_heap_permutation(const PathHandleGraph& graph,
                               const std::vector<std::vector<path_handle_t>>& path_groups,
                               uint64_t n_permutations,
                               const std::function<void(const std::vector<uint64_t>&, uint64_t)>& func);

}

}
