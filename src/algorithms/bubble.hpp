#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <set>
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
#include <atomic_bitvector.hpp>

namespace odgi {

using namespace handlegraph;

namespace algorithms {

void for_each_bubble(const PathHandleGraph& graph,
                     const path_handle_t& path,
                     const std::function<void(const step_handle_t& begin,
                                              const step_handle_t& end)>& func);

}

}
