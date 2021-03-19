#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "hash_map.hpp"
#include "odgi.hpp"

namespace odgi {
    namespace algorithms {

        // BED files are 0-based http://genome.ucsc.edu/FAQ/FAQformat#format1
        struct min_bed_record_t {
            uint64_t chromStart;
            uint64_t chromEnd;
            double path_layout_nuc_dist_ratio;
        };

        ska::flat_hash_map<path_handle_t, std::vector<min_bed_record_t>> find_ranges(const std::string &bed_in_file,
                                                                               const double &min_median_factor = 3,
                                                                               const uint64_t &relax_num = 10,
                                                                               const double &relax_percentage = 0.1,
                                                                               const bool &relax_num_ = true);

        min_bed_record_t parse_bed_line(const std::string &bed_line, const char delim);

    }

}