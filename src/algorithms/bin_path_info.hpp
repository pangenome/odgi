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

struct pathinfo_t {
    double mean_cov;
    double mean_inv;
    double mean_pos;
};

void bin_path_info(const PathHandleGraph& graph,
                   const std::string& prefix_delimiter,
                   const uint64_t& num_bins,
                   std::vector<std::pair<std::string, std::vector<pathinfo_t>>>& table);

}
}
