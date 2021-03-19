#include "relax_prep_bed_records.hpp"

namespace odgi {
    namespace algorithms {

        ska::flat_hash_map<path_handle_t, std::vector<min_bed_record_t>> find_ranges(const std::string &bed_in_file,
                                                                               const double &min_median_factor,
                                                                               const uint64_t &relax_num,
                                                                               const double &relax_percentage,
                                                                               const bool &relax_num_b) {
            std::ifstream bed_in_stream(bed_in_file);
            std::string bed_line;
            bool first_line_was_read = false;
            const char delim = '\t';
            ska::flat_hash_map<path_handle_t, vector<min_bed_record_t>> min_bed_record_t_map;
            std::vector<min_bed_record_t> min_bed_record_t_vec;
            while (std::getline(bed_in_stream, bed_line)) {
                if (!first_line_was_read) {
                    first_line_was_read = true;
                }
                min_bed_record_t_vec.push_back(parse_bed_line(bed_line, delim));
            }
            for (auto &min_bed_record : min_bed_record_t_vec) {
                std::cout << min_bed_record.chromStart << "\t" << min_bed_record.chromEnd << "\t" << min_bed_record.path_layout_nuc_dist_ratio << std::endl;
            }
            return min_bed_record_t_map;
        }

        min_bed_record_t parse_bed_line(const std::string &bed_line, const char delim) {
            uint8_t tab_pos = 0;
            std::string token;
            std::istringstream tokenStream(bed_line);
            min_bed_record_t m_b_r_t;
            while (std::getline(tokenStream, token, delim)) {
                 if (tab_pos == 1) {
                    m_b_r_t.chromStart = std::stoi(token);
                } else if (tab_pos == 2) {
                    m_b_r_t.chromEnd = std::stoi(token);
                } else if (tab_pos == 5) {
                    m_b_r_t.path_layout_nuc_dist_ratio = std::stod(token);
                }
                tab_pos++;
            }
            return m_b_r_t;
        }

    }

}