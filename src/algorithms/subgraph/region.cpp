#include <iostream>
#include <fstream>
#include <cassert>
#include <odgi.hpp>
#include <position.hpp>
#include "region.hpp"
#include "split.hpp"

namespace odgi {

    void parse_region(const std::string& target, std::string& name, int64_t& start, int64_t& end) {
        start = -1;
        end = -1;
        const size_t foundLastColon = target.find_last_of(":");

        // we only have a single string, use the whole sequence as the target
        if (foundLastColon == std::string::npos) {
            name = target;
        } else {
            name = target.substr(0, foundLastColon);
            size_t foundRangeDash = target.find("-", foundLastColon);
            if (foundRangeDash == std::string::npos) {
                start = atoi(target.substr(foundLastColon + 1).c_str());
                end = start;
            } else {
                start = atoi(target.substr(foundLastColon + 1, foundRangeDash - foundRangeDash - 1).c_str());
                end = atoi(target.substr(foundRangeDash + 1).c_str());
            }
        }
    }

    void parse_bed_regions(const std::string& bed_path,
                           std::vector<Region>& out_regions,
                           std::vector<std::string>* out_names) {
        out_regions.clear();
        std::ifstream bedstream(bed_path);
        if (!bedstream) {
            std::cerr << "[odgi::parse_bed_regions] unable to open bed file: " << bed_path << std::endl;
            return;
        }
        std::string row;
        std::string sbuf;
        std::string ebuf;
        std::string nbuf;
        for (int line = 1; getline(bedstream, row); ++line) {
            Region region;
            if (row.size() < 2 || row[0] == '#') {
                continue;
            }
            std::istringstream ss(row);
            if (!getline(ss, region.seq, '\t') ||
                !getline(ss, sbuf, '\t') ||
                !getline(ss, ebuf, '\t') ||
                (out_names != nullptr && !getline(ss, nbuf, '\t'))) {
                std::cerr << "[odgi::parse_bed_regions] error parsing bed line " << line << ": " << row << std::endl;
            } else {
                region.start = std::stoi(sbuf);
                region.end = std::stoi(ebuf);
                assert(region.end > region.start);

                // convert from BED-style to 0-based inclusive coordinates
                region.end -= 1;

                out_regions.push_back(region);

                if (out_names != nullptr) {
                    out_names->push_back(nbuf);
                }
            }
        }
    }

    void add_bed_range(std::vector<odgi::path_range_t>& path_ranges,
                       const odgi::graph_t &graph,
                       const std::string &buffer) {
        if (!buffer.empty() && buffer[0] != '#') {
            const auto vals = split(buffer, '\t');

            const auto &path_name = vals[0];
            if (!graph.has_path(path_name)) {
                std::cerr << "[odgi::add_bed_range] error: path " << path_name << " not found in graph" << std::endl;
                exit(1);
            } else {
                const uint64_t start = vals.size() > 1 ? (uint64_t) std::stoi(vals[1]) : 0;
                uint64_t end = 0;
                if (vals.size() > 2) {
                    end = (uint64_t) std::stoi(vals[2]);
                } else {
                    // In the BED format, the end is non-inclusive, unlike start
                    graph.for_each_step_in_path(graph.get_path_handle(path_name), [&](const step_handle_t &s) {
                        end += graph.get_length(graph.get_handle_of_step(s));
                    });
                }

                if (start >= end) {
                    std::cerr << "[odgi::add_bed_range] error: wrong input coordinates in row: " << buffer << std::endl;
                    exit(1);
                }

                path_ranges.push_back(
                        {
                            {
                                graph.get_path_handle(path_name),
                                start,
                                false
                            },
                            {
                                graph.get_path_handle(path_name),
                                end,
                                false
                            },
                            (vals.size() > 5 && vals[5] == "-"),
                            vals.size() > 3 ? vals[3] : ".",
                            buffer
                        });
            }
        }
    };

}