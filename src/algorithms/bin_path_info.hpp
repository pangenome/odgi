#pragma once

#include <vector>
#include <utility>
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/path_handle_graph.hpp>

namespace odgi {
namespace algorithms {

using namespace std;
using namespace handlegraph;

void bin_path_coverage(const PathHandleGraph& graph,
                       const std::string& prefix_delimiter,
                       const uint64_t& num_bins,
                       std::vector<std::pair<std::string, std::vector<double>>>& table);

}
}
