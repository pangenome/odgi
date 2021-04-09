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
            std::string cur_chrom;
            while (std::getline(bed_in_stream, bed_line)) {
                if (!first_line_was_read) {
                    first_line_was_read = true;
                    vector<basic_string<char>> vals = split(bed_line, '\t');
                    cur_chrom = vals[0];
                }
                min_bed_record_t_vec.push_back(parse_bed_line(bed_line, delim)); // TODO parse_bed_line must get the splitted bed_line as argument so we can figure out if the CHROM changed!
                // TODO TAKE OVER DESCRIPTION FROM relax_main

            }
            for (auto &min_bed_record : min_bed_record_t_vec) {
                std::cout << min_bed_record.chromStart << "\t" << min_bed_record.chromEnd << "\t" << min_bed_record.path_layout_nuc_dist_ratio << std::endl;
            }
            return min_bed_record_t_map;
        }

        min_bed_record_t parse_bed_line(const std::string &bed_line, const char delim) {
            auto vals = split(bed_line, '\t');
            uint8_t tab_pos = 0;
            min_bed_record_t m_b_r_t;
            for (auto &val : vals) {
                if (tab_pos == 1) {
                    m_b_r_t.chromStart = std::stoi(val);
                } else if (tab_pos == 2) {
                    m_b_r_t.chromEnd = std::stoi(val);
                } else if (tab_pos == 5) {
                    m_b_r_t.path_layout_nuc_dist_ratio = std::stod(val);
                }
                tab_pos++;
            }
            return m_b_r_t;
        }

    }

}