#pragma once

#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <vector>
#include <set>
#include <deque>
#include "hash_map.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

void diff_priv(
    const PathHandleGraph& input,
    PathHandleGraph& priv,
    double target_coverage,
    uint64_t bp_limit = 0);

}
}
