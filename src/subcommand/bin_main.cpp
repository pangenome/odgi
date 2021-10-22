#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/bin_path_info.hpp"
#include "algorithms/bin_path_depth.hpp"
#include "gfa_to_handle.hpp"
#include "utils.hpp"

#include <regex>

namespace odgi {

using namespace odgi::subcommand;

int main_bin(int argc, char** argv) {

    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi bin";
    argv[0] = (char*)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("Binning of pangenome sequence and path information in the graph.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this FILE. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::Group bin_opts(parser, "[ Bin Options ]");
    args::ValueFlag<std::string> path_delim(bin_opts, "path-delim", "Annotate rows by prefix and suffix of this delimiter.", {'D', "path-delim"});
    args::Flag aggregate_delim(bin_opts, "aggregate-delim", "Aggregate on path prefix delimiter. Argument depends on -D,--path-delim=[STRING].", {'a', "aggregate-delim"});
    args::Flag output_json(bin_opts, "write-json", "Print bins and links to stdout in pseudo JSON format. "
                                                 "Each line is a  valid JSON object, but the whole file is not a valid JSON! "
                                                 "First, each  bin including its pangenome sequence is printed to stdout per line.  "
                                                 "Second, for each path in the graph, its traversed bins including  metainformation: **bin** (bin identifier), "
                                                 "**mean.cov** (mean coverage  of the path in this bin), "
                                                 "**mean.inv** (mean inversion rate of this  path in this bin), "
                                                 "**mean.pos** (mean nucleotide position of this path  in this bin), "
                                                 "and an array of ranges determining the nucleotide  position of the path in this bin. "
                                                 "Switching first and last nucleotide  in a range represents a complement "
                                                 "reverse orientation of that  particular sequence.", {'j', "json"});
    args::ValueFlag<uint64_t> num_bins(bin_opts, "N", "The number of bins the pangenome sequence should be chopped up to.", {'n', "num-bins"});
    args::ValueFlag<uint64_t> bin_width(bin_opts, "bp", "The bin width specifies the size of each bin.", {'w', "bin-width"});
    args::Flag write_seqs_not(bin_opts, "write-seqs-not", "If -j,--json is set, no nucleotide sequences will be printed to stdout in order to save disk space.", {'s', "no-seqs"});
    args::Flag drop_gap_links(bin_opts, "drop-gap-links", "Don't include gap links in the output. "
                                                          "We divide links into 2 classes:\n1. The links which help to follow complex variations. "
                                                          "They need to be drawn, else one could not follow the sequence of a path.\n"
                                                          "2. The links helping to follow simple variations. These links are called gap-links."
                                                          " Such links solely connecting a path from left to right may not be"
                                                          "relevant to understand a path's traveral through the bins. Therfore,"
                                                          " when this option is set, the gap-links are left out saving disk space.", {'g', "no-gap-links"});
    args::Group haplo_blocker_opts(parser, "[ HaploBlocker Options ]");
    args::Flag haplo_blocker(haplo_blocker_opts, "haplo-blocker", "Write a TSV to stdout formatted in a "
                                                                  "way ready for HaploBlocker: Each row corresponds to a node. "
                                                                  "Each column corresponds to a path. Each value is the coverage of "
                                                                  "a specific node of a specific path.", {'b', "haplo-blocker"});
    args::ValueFlag<uint64_t> haplo_blocker_min_paths(haplo_blocker_opts, "N", "Specify the minimum number of paths "
                                                                               "that need to be present in the bin to actually"
                                                                               " report that bin (default: 1).", {'p', "haplo-blocker-min-paths"});
    args::ValueFlag<uint64_t> haplo_blocker_min_depth(haplo_blocker_opts, "N", "Specify the minimum depth a path needs to have in a bin to actually report that bin (default: 1).", {'c', "haplo-blocker-min-depth"});
	args::Group threading(parser, "[ Threading ]");
	args::ValueFlag<uint64_t> nthreads(threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
    args::Group processing_info_opts(parser, "[ Processing Information ]");
    args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi bin.", {'h', "help"});
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    if (!dg_in_file) {
        std::cerr << "[odgi::bin] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

	const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

	graph_t graph;
    assert(argc > 0);
    if (!args::get(dg_in_file).empty()) {
        std::string infile = args::get(dg_in_file);
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
			utils::handle_gfa_odgi_input(infile, "bin", args::get(progress), num_threads, graph);
        }
    }

    std::string delim = args::get(path_delim);
    bool agg_delim = args::get(aggregate_delim);
    auto get_path_prefix = [&](const std::string& path_name) -> std::string {
        if (agg_delim || delim.empty()) {
            return "NA";
        } else {
            return path_name.substr(0, path_name.find(delim));
        }
    };
    auto get_path_suffix = [&](const std::string& path_name) -> std::string {
        if (agg_delim || delim.empty()) {
            return "NA";
        } else {
            return path_name.substr(path_name.find(delim)+1);
        }
    };

    // our aggregation matrix
    std::vector<std::pair<std::string, std::vector<algorithms::path_info_t>>> table;
    if (args::get(num_bins) + args::get(bin_width) == 0) {
        std::cerr << "[odgi::bin] error: a bin width or a bin count is required" << std::endl;
        return 1;
    }

    if (haplo_blocker) {
        std::cerr << "[odgi::bin] main: running in HaploBlocker mode. Ignoring input parameters -D/--path-delim, -j/--json, -a/--aggregate-delim, "
                     "-n/--num-bins, -w/--bin-width, -s/--no-seqs, -g/--no-gap-links." << std::endl;
        // first pass: collect #nucleotides, fill in_all_bins_bv, unique_bins_bv
        uint64_t haplo_blocker_min_paths_ = args::get(haplo_blocker_min_paths) ? args::get(haplo_blocker_min_paths) : 1;
        uint64_t haplo_blocker_min_depth_ = args::get(haplo_blocker_min_depth) ? args::get(haplo_blocker_min_depth) : 1.0;
        algorithms::bin_path_depth(graph, args::get(progress),
                                      haplo_blocker_min_paths_, haplo_blocker_min_depth_);
        // write header of table to stdout

        // for all paths, for each node and nucleotide in path --> better unordered_map https://stackoverflow.com/questions/1939953/how-to-find-if-a-given-key-exists-in-a-c-stdmap
        // map<uint64_t, double>: ++ when we cover the bin
        // as it is done so far

        // update mean depth
        // write to std::cout, only if we have !in_all_bins_bv[idx] && !unique_bins_bv[idx] && auto it = m.find("f"); if (it != m.end()) {/*Use it->second*/}.
        // cap depth by 255
    } else {

        // ODGI JSON VERSION
        const uint64_t ODGI_JSON_VERSION = 12; // this brings the exact nucleotide positions for each bin for each path referred to as ranges

        std::function<void(const uint64_t&, const uint64_t&)> write_header_tsv
                = [&] (const uint64_t pangenome_length, const uint64_t bin_width) {
                    // no header necessary for tsv so far
                };

        std::function<void(const uint64_t&,
                           const uint64_t&)> write_header_json
                = [&] (const uint64_t pangenome_length, const uint64_t bin_width) {
                    std::cout << "{\"odgi_version\": " << ODGI_JSON_VERSION << ",";
                    std::cout << "\"bin_width\": " << bin_width << ",";
                    std::cout << "\"pangenome_length\": " << pangenome_length << "}" << std::endl;
                };

        std::function<void(const uint64_t&,
                           const std::string&)> write_seq_json
                = [&](const uint64_t& bin_id, const std::string& seq) {
                    if (args::get(write_seqs_not)) {
                        std::cout << "{\"bin_id\":" << bin_id << "}" << std::endl;
                    } else {
                        std::cout << "{\"bin_id\":" << bin_id << ","
                                  << "\"sequence\":\"" << seq << "\"}" << std::endl;
                    }
                };

        std::function<void(const vector<std::pair<uint64_t , uint64_t >>&)> write_ranges_json
                = [&](const vector<std::pair<uint64_t , uint64_t >>& ranges) {
                    std::cout << "[";
                    for (int i = 0; i < ranges.size(); i++) {
                        std::pair<uint64_t, uint64_t > range = ranges[i];
                        if (i == 0) {
                            std::cout << "[" << range.first << "," << range.second << "]";
                        } else {
                            std::cout << "," << "[" << range.first << "," << range.second << "]";
                        }
                    }
                    std::cout << "]";
                };

        std::function<void(const std::string&,
                           const std::vector<std::pair<uint64_t, uint64_t>>&,
                           const std::map<uint64_t, algorithms::path_info_t>&)> write_json
                = [&](const std::string& path_name,
                      const std::vector<std::pair<uint64_t, uint64_t>>& links,
                      const std::map<uint64_t, algorithms::path_info_t>& bins) {
                    std::string name_prefix = get_path_prefix(path_name);
                    std::string name_suffix = get_path_suffix(path_name);
                    std::cout << R"({"path_name":")" << path_name << "\",";
                    if (!delim.empty()) {
                        std::cout << "\"path_name_prefix\":\"" << name_prefix << "\","
                                  << "\"path_name_suffix\":\"" << name_suffix << "\",";
                    }
                    std::cout << "\"bins\":[";
                    auto entry_it = bins.begin();
                    for (uint64_t i = 0; i < bins.size(); ++i) {
                        auto& bin_id = entry_it->first;
                        auto &info = entry_it->second;
                        std::cout << "[" << bin_id << ","
                                  << info.mean_depth << ","
                                  << info.mean_inv << ","
                                  << info.mean_pos << ",";
                        write_ranges_json(info.ranges);
                        std::cout << "]";
                        if (i+1 != bins.size()) {
                            std::cout << ",";
                        }
                        ++entry_it;
                    }
                    std::cout << "]";
                    std::cout << ",\"links\":[";
                    for (uint64_t i = 0; i < links.size(); ++i) {
                        auto &link = links[i];
                        std::cout << "[" << link.first << "," << link.second << "]";
                        if (i + 1 < links.size()) std::cout << ",";
                    }
                    std::cout << "]}" << std::endl;
                };

        std::function<void(const uint64_t&,
                           const std::string&)> write_seq_noop
                = [&](const uint64_t& bin_id, const std::string& seq) {
                };

        std::function<void(const std::string&,
                           const std::vector<std::pair<uint64_t, uint64_t>>&,
                           const std::map<uint64_t, algorithms::path_info_t>&)> write_tsv
                = [&](const std::string& path_name,
                      const std::vector<std::pair<uint64_t, uint64_t>>& links,
                      const std::map<uint64_t, algorithms::path_info_t>& bins) {
                    std::string name_prefix = get_path_prefix(path_name);
                    std::string name_suffix = get_path_suffix(path_name);
                    for (auto& entry : bins) {
                        auto& bin_id = entry.first;
                        auto& info = entry.second;
                        if (info.mean_depth > 0) {
                            std::cout << path_name << "\t"
                                      << name_prefix << "\t"
                                      << name_suffix << "\t"
                                      << bin_id << "\t"
                                      << info.mean_depth << "\t"
                                      << info.mean_inv << "\t"
                                      << info.mean_pos << "\t"
                                      << info.ranges[0].first << "\t";
                            if (info.ranges[info.ranges.size() - 1].second == 0) {
                                std::cout << info.ranges[info.ranges.size() - 1].first << std::endl;
                            } else {
                                std::cout << info.ranges[info.ranges.size() - 1].second << std::endl;
                            }
                        }
                    }
                };

        if (args::get(output_json)) {
            algorithms::bin_path_info(graph, (args::get(aggregate_delim) ? args::get(path_delim) : ""),
                                      write_header_json,write_json, write_seq_json,
                                      args::get(num_bins), args::get(bin_width), args::get(drop_gap_links),
                                      args::get(progress));
        } else {
            std::cout << "path.name" << "\t"
                      << "path.prefix" << "\t"
                      << "path.suffix" << "\t"
                      << "bin" << "\t"
                      << "mean.cov" << "\t"
                      << "mean.inv" << "\t"
                      << "mean.pos" << "\t"
                      << "first.nucl" << "\t"
                      << "last.nucl" << std::endl;
            algorithms::bin_path_info(graph, (args::get(aggregate_delim) ? args::get(path_delim) : ""),
                                      write_header_tsv,write_tsv, write_seq_noop,
                                      args::get(num_bins), args::get(bin_width), args::get(drop_gap_links),
                                      args::get(progress));
        }
    }
    return 0;
}

static Subcommand odgi_bin("bin", "Binning of pangenome sequence and path information in the graph.",
                              PIPELINE, 3, main_bin);


}
