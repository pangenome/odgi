#pragma once

#include <vector>
#include <utility>
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <iomanip> // std::setprecision
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include "progress.hpp"

namespace odgi {
    namespace algorithms {

        using namespace std;
        using namespace handlegraph;

        struct path_info_t {
            double mean_depth;
            double mean_inv;
            double mean_pos;
            vector<std::pair<uint64_t,uint64_t>> ranges;
            // long int first_nucleotide;
            // long int last_nucleotide;
        };

        void bin_path_info(const PathHandleGraph &graph,
                           const std::string &prefix_delimiter,
                           const std::function<void(const uint64_t &, const uint64_t &)> &handle_header,
                           const std::function<void(const std::string &,
                                                    const std::vector<std::pair<uint64_t, uint64_t>> &,
                                                    const std::map<uint64_t, algorithms::path_info_t> &)> &handle_path,
                           const std::function<void(const uint64_t &, const std::string &)> &handle_sequence,
                           uint64_t num_bins = 0,
                           uint64_t bin_width = 0,
                           bool drop_gap_links = false,
                           bool progress = false);
    }
}
