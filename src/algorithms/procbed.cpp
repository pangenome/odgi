#include "procbed.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

void adjust_ranges(const PathHandleGraph& graph, const std::string& bed_targets) {

    // collect the subgraph path map
    ska::flat_hash_map<std::string, std::vector<interval_t>> subpaths;
    // iterate over paths
    graph.for_each_path_handle(
        [&subpaths,&graph](const path_handle_t& path) {
            // check if the path is named following pansn
            auto name = graph.get_path_name(path);
            std::string base;
            uint64_t start = 0;
            uint64_t end = 0;
            auto c = name.find(':');
            auto d = name.find('-', c);
            if (c != std::string::npos && d != std::string::npos) {
                // PanSN
                // if so, collect its name and length and try to put it into our subpath
                base = name.substr(0,c);
                start = std::stoul(name.substr(c+1,d));
                end = std::stoul(name.substr(d+1));
            } else {
                // if not, measure its length and use [0, length) as our interval
                base = name;
                uint64_t len = 0;
                graph.for_each_step_in_path(
                    path,
                    [&graph,&len](const step_handle_t& step) {
                        len += graph.get_length(graph.get_handle_of_step(step));
                    });
                start = 0;
                end = len;
            }
            subpaths[base].push_back(interval_t(start, end));
        });
    // sort the intervals
    for (auto& p : subpaths) {
        std::sort(p.second.begin(), p.second.end());
    }

    ska::flat_hash_map<std::string, std::vector<std::pair<interval_t, std::string>>> bed_intervals;
    std::ifstream bed(bed_targets.c_str());
    std::string line;
    while (std::getline(bed, line)) {
        if (!line.empty()) {
            auto vals = split(line, '\t');
            if (vals.size() < 4) {
                std::cerr << "[odgi::algorithms::adjust_ranges]"
                          << "BED line does not have enough fields to define an interval"
                          << std::endl << line << std::endl;
                std::abort();
            }
            auto& path_name = vals[0];
            uint64_t start = std::stoul(vals[1]);
            uint64_t end = std::stoul(vals[2]);
            const std::string& name = vals[3];
            bed_intervals[path_name].push_back(make_pair(interval_t(start, end), name));
        }
    }
    // sort the intervals
    for (auto& p : bed_intervals) {
        std::sort(p.second.begin(), p.second.end());
    }

    // now we match BED ranges to graph intervals
    // using a two-list sweep to find ranges that fit into our graph
    // we can either emit warnings for those that are intersected
    // or we can cut them
    for (auto& b : bed_intervals) {
        auto& ref = b.first;
        auto& bedivals = b.second;
        auto f = subpaths.find(ref);
        if (f != subpaths.end()) {
            auto& refivals = f->second;
            // map from range end to starting positions
            std::map<uint64_t, std::vector<uint64_t>> ref_range_ends;
            auto r = refivals.begin();
            auto b = bedivals.begin();
            while (b != bedivals.end() && (r != refivals.end() || !ref_range_ends.empty())) {
                auto& b_start = b->first.first;
                auto& b_end = b->first.second;
                auto& b_key = b->second;
                while (ref_range_ends.size()
                       && ref_range_ends.begin()->first < b_start) {
                    ref_range_ends.erase(ref_range_ends.begin());
                }
                while (r != refivals.end() && r->second < b_start) {
                    // non-overlapping
                    ++r;
                }
                while (r != refivals.end() && r->first < b_end) {
                    ref_range_ends[r->second].push_back(r->first);
                    ++r;
                }
                // now we check if b is in the open ranges
                // and for each one we'll do a mapping
                for (auto& f : ref_range_ends) {
                    auto& ref_end = f.first;
                    if (ref_end >= b_end) {
                        // find the ranges that can contain this interval
                        for (auto& ref_start : f.second) {
                            if (b_start >= ref_start && b_end > ref_start) {
                                std::cout << ref << ":" << ref_start << "-" << ref_end << "\t"
                                          << b_start - ref_start << "\t" << b_end - ref_start << "\t"
                                          << b_key << std::endl;
                            }
                        }
                    }
                }
                ++b;
            }
        }
    }
}

}

}
