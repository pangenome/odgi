#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/bin_path_info.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_bin(int argc, char** argv) {

    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi bin";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("binning of path information in the graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the graph in this file", {'o', "out"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> path_delim(parser, "path-delim", "annotate rows by prefix and suffix of this delimiter", {'D', "path-delim"});
    args::Flag output_json(parser, "write-json", "write JSON format output including additional path positional information", {'j', "json"});
    args::Flag aggregate_delim(parser, "aggregate-delim", "aggregate on path prefix delimiter", {'a', "aggregate-delim"});
    args::ValueFlag<uint64_t> num_bins(parser, "N", "number of bins", {'n', "num-bins"});
    args::ValueFlag<uint64_t> bin_width(parser, "bp", "width of each bin in basepairs along the graph vector", {'w', "bin-width"});
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

    graph_t graph;
    assert(argc > 0);
    std::string infile = args::get(dg_in_file);
    if (infile.size()) {
        if (infile == "-") {
            graph.load(std::cin);
        } else {
            ifstream f(infile.c_str());
            graph.load(f);
            f.close();
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
        std::cerr << "[odgi bin] error: a bin width or a bin count is required" << std::endl;
        return 1;
    }

    std::function<void(const std::string&, const hash_map<uint64_t, algorithms::path_info_t>&)> write_json
        = [&](const std::string& path_name,
              const hash_map<uint64_t, algorithms::path_info_t>& table) {
        std::string name_prefix = get_path_prefix(path_name);
        std::string name_suffix = get_path_suffix(path_name);
        for (auto& entry : table) {
            auto& bin_id = entry.first;
            auto& info = entry.second;
            if (info.mean_cov) {
                std::cout << "{\"path_name\":\"" << path_name << "\",";
                if (!delim.empty()) {
                    std::cout << "\"path_name_prefix\":\"" << name_prefix << "\","
                              << "\"path_name_suffix\":\"" << name_suffix << "\",";
                }
                std::cout << "\"bin_id\":" << bin_id << ","
                          << "\"mean_cov\":" << info.mean_cov << ","
                          << "\"mean_inv\":" << info.mean_inv << ","
                          << "\"mean_pos\":" << info.mean_pos << ","
                          << "\"begins\":[";
                for (uint64_t i = 0; i < info.begins.size(); ++i) {
                    auto& p = info.begins[i];
                    std::cout << "[" << p.other_bin << ","
                              << p.pos_in_other_bin << ","
                              << p.pos_in_this_bin << ","
                              << p.pos_in_path << ","
                              << (p.is_rev?"true":"false") << "]";
                    if (i+1 < info.begins.size()) {
                        std::cout << ",";
                    }
                }
                std::cout << "],"
                          << "\"ends\":[";
                for (uint64_t i = 0; i < info.ends.size(); ++i) {
                    auto& p = info.ends[i];
                    std::cout << "[" << p.other_bin << ","
                              << p.pos_in_other_bin << ","
                              << p.pos_in_this_bin << ","
                              << p.pos_in_path << ","
                              << (p.is_rev?"true":"false") << "]";
                    if (i+1 < info.ends.size()) {
                        std::cout << ",";
                    }
                }
                std::cout << "]}" << std::endl;
            }
        }
    };

    std::function<void(const std::string&, const hash_map<uint64_t, algorithms::path_info_t>&)> write_tsv
        = [&](const std::string& path_name,
              const hash_map<uint64_t, algorithms::path_info_t>& table) {
        std::string name_prefix = get_path_prefix(path_name);
        std::string name_suffix = get_path_suffix(path_name);
        for (auto& entry : table) {
            auto& bin_id = entry.first;
            auto& info = entry.second;
            if (info.mean_cov) {
                std::cout << path_name << "\t"
                          << name_prefix << "\t"
                          << name_suffix << "\t"
                          << bin_id << "\t"
                          << info.mean_cov << "\t"
                          << info.mean_inv << "\t"
                          << info.mean_pos << "\t" << std::endl;
            }
        }
    };

    if (args::get(output_json)) {
        algorithms::bin_path_info(graph, (args::get(aggregate_delim) ? args::get(path_delim) : ""),
                                  write_json,
                                  args::get(num_bins), args::get(bin_width));
    } else {
        std::cout << "path.name" << "\t" << "path.prefix" << "\t" << "path.suffix" << "\t" << "bin" << "\t" << "mean.cov" << "\t" << "mean.inv" << "\t" << "mean.pos" << std::endl;
        algorithms::bin_path_info(graph, (args::get(aggregate_delim) ? args::get(path_delim) : ""),
                                  write_tsv,
                                  args::get(num_bins), args::get(bin_width));
    }

    return 0;
}

static Subcommand odgi_bin("bin", "bin path information across the graph",
                              PIPELINE, 3, main_bin);


}
