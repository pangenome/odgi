#pragma once

#include <vector>
#include <utility>
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/path_handle_graph.hpp>

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

        typedef std::pair<uint64_t, uint64_t> link_t;
        typedef std::vector<link_t> link_vec_t;
        typedef std::map<uint64_t, algorithms::path_info_t> bin_map_t;

        struct BinSerializer {
            typedef algorithms::link_vec_t link_vec_t;
            typedef algorithms::bin_map_t bin_map_t;

            std::string path_delim;
            bool aggregate_delim;

            std::string get_path_prefix(const std::string& path_name);
            std::string get_path_suffix(const std::string& path_name);

            BinSerializer(const std::string& path_delim, bool aggregate_delim);
            virtual ~BinSerializer();

            virtual void write_header(const uint64_t pangenome_length,
                                      const uint64_t bin_width) = 0;

            virtual void write_seq(const uint64_t& bin_id, const std::string& seq) = 0;

            virtual void write_path(const std::string& path_name,
                                    const link_vec_t&, const bin_map_t&) = 0;
        };

        void bin_path_info(const PathHandleGraph &graph,
                           const std::string &prefix_delimiter,
                           std::shared_ptr<BinSerializer>& bin_serializer,
                           const std::function<void(const std::string&)> &handle_fasta,
                           uint64_t num_bins = 0,
                           uint64_t bin_width = 0,
                           bool drop_gap_links = false);
    }
}
