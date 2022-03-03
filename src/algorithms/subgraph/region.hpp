#ifndef ODGI_REGION_H
#define ODGI_REGION_H

#include <string>
#include <vector>
#include <sstream>
#include "position.hpp"

namespace odgi {

    // Represent a parsed genomic region.
    // A -1 for start or end indicates that that coordinate is not used.
    // Generally regions parsed form user input will be 1-based.
    struct Region {
        std::string seq;
        int64_t start;
        int64_t end;
    };

    // Parse a genomic contig[:start-end] region. Outputs -1 for missing start or end.
    void parse_region(const std::string &target, std::string &name, int64_t &start, int64_t &end);


    // Parse a genomic contig[:start-end] region. Outputs -1 for missing start or end.
    inline void parse_region(std::string &region, Region &out_region) {
        parse_region(region,
                     out_region.seq,
                     out_region.start,
                     out_region.end);
    }

    // parse a bed file and return a list of regions (like above)
    // IMPORTANT: expects usual 0-based BED format.
    // So bedline "chr1   5   10" will return start=5 stop=9
    void parse_bed_regions(
            const std::string &bed_path,
            std::vector<Region> &out_regions,
            std::vector<std::string> *out_names = nullptr);

    void add_bed_range(std::vector<odgi::path_range_t>& path_ranges,
                       const odgi::graph_t &graph,
                       const std::string &buffer);
}

#endif //ODGI_REGION_H
