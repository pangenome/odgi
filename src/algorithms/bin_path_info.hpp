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
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/path_handle_graph.hpp>

#include "sdsl/coder_elias_delta.hpp"
#include "sdsl/coder_elias_gamma.hpp"
#include "sdsl/coder_fibonacci.hpp"
#include "sdsl/dac_vector.hpp"
#include "sdsl/enc_vector.hpp"
#include "sdsl/vlc_vector.hpp"

namespace odgi {
namespace algorithms {

using namespace std;
using namespace handlegraph;

struct path_info_t {
    double mean_cov;
    double mean_inv;
    double mean_pos;
    long int first_nucleotide;
    long int last_nucleotide;
};

void bin_path_info(const PathHandleGraph& graph,
                   const std::string& prefix_delimiter,
                   const std::function<void(const std::string&,
                                            const std::vector<std::pair<uint64_t, uint64_t>>&,
                                            const std::map<uint64_t, algorithms::path_info_t>&)>& handle_path,
                   const std::function<void(const uint64_t&, const std::string&)>& handle_sequence,
                   // const std::function<void(const std::map<std::string&, std::vector<uint64_t&>>)> handle_index, // TODO @ekg Is this correct like that?
                   uint64_t num_bins = 0,
                   uint64_t bin_width = 0,
                   bool index_json = false,
                   const std::string& eval_out_file = "");
}
}
