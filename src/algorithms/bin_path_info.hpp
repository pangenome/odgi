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
#include "../hash_map.hpp"

namespace odgi {
namespace algorithms {

using namespace std;
using namespace handlegraph;

struct path_info_t {
    struct path_pos_t {
        int64_t other_bin;
        uint64_t pos_in_other_bin;
        uint64_t pos_in_this_bin;
        uint64_t pos_in_path;
        bool is_rev;
    };
    double mean_cov;
    double mean_inv;
    double mean_pos;
    std::vector<path_pos_t> begins;
    std::vector<path_pos_t> ends;
};

void bin_path_info(const PathHandleGraph& graph,
                   const std::string& prefix_delimiter,
                   const std::function<void(const std::string&, const hash_map<uint64_t, algorithms::path_info_t>&)>& handle_path,
                   uint64_t num_bins = 0,
                   uint64_t bin_width = 0);


}
}
