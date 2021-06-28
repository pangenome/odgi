#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <omp.h>
#include <handlegraph/types.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include "flat_hash_map.hpp"

#include "perfect_neighbors.hpp"
#include "dset64.hpp"
#include "BooPHF.h"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

typedef boomphf::mphf<uint64_t, boomphf::SingleHashFunctor<uint64_t>> boophf_uint64_t;

std::vector<std::vector<handle_t>> simple_components(
    const PathHandleGraph &graph, const uint64_t& min_size, const bool& return_all_handles, const uint64_t& nthreads);

}

}
