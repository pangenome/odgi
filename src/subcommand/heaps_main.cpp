#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/heaps.hpp"
#include "utils.hpp"
#include "split.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_heaps(int argc, char **argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc - 1; ++i) {
        argv[i] = argv[i + 1];
    }
    const std::string prog_name = "odgi heaps";
    argv[0] = (char *) prog_name.c_str();
    --argc
;
    args::ArgumentParser parser("Extract matrix of path pangenome coverage permutations for power law regression.");
    args::Group mandatory_opts(parser, "[ MANDATORY ARGUMENTS ]");
    args::ValueFlag<std::string> og_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::Group heaps_opts(parser, "[ Heaps Options ]");
    args::ValueFlag<std::string> _path_groups(heaps_opts, "FILE", "Group paths as described in two-column FILE, with columns path.name and group.name.",
                                              {'p', "path-groups"});
    args::Flag group_by_sample(heaps_opts, "bool", "Following PanSN naming (sample#hap#ctg), group by sample (1st field).", {'S', "group-by-sample"});
    args::Flag group_by_haplotype(heaps_opts, "bool", "Following PanSN naming (sample#hap#ctg), group by haplotype (2nd field).", {'H', "group-by-haplotype"});
    args::ValueFlag<std::string> _bed_targets(heaps_opts, "FILE", "BED file over path space of the graph, describing a subset of the graph to consider.",
                                              {'b', "bed-targets"});
    args::ValueFlag<uint64_t> _n_permutations(heaps_opts, "N", "Number of permutations to run.",
                                             {'n', "n-permutations"});
    args::ValueFlag<uint64_t> _min_node_depth(heaps_opts, "N", "Exclude nodes with less than this path depth (default: 0).",
                                         {'d', "min-node-depth"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> nthreads(threading_opts, "N", "Number of threads to use for parallel operations.",
                                       {'t', "threads"});
    args::Group processing_info_opts(parser, "[ Processing Information ]");
    //args::Flag debug(processing_info_opts, "debug", "Print information about the process to stderr.", {'d', "debug"});
    args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi heaps.", {'h', "help"});
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
    if (argc == 1) {
        std::cout << parser;
        return 1;
    }

    if (!og_in_file) {
        std::cerr
            << "[odgi::heaps] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
            << std::endl;
        return 1;
    }

    const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;
    omp_set_num_threads(num_threads);

    const uint64_t n_permutations = args::get(_n_permutations) ? args::get(_n_permutations) : 1;

    const uint64_t min_node_depth = args::get(_min_node_depth) ? args::get(_min_node_depth) : 0;

    graph_t graph;
    assert(argc > 0);
    {
        const std::string infile = args::get(og_in_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
                utils::handle_gfa_odgi_input(infile, "heaps", args::get(progress), num_threads, graph);
            }
        }
    }

    std::vector<std::vector<path_handle_t>> path_groups;
    {
        ska::flat_hash_map<std::string, std::vector<path_handle_t>> path_groups_map;
        if (_path_groups) {
            std::ifstream refs(args::get(_path_groups).c_str());
            std::string line;
            while (std::getline(refs, line)) {
                if (!line.empty()) {
                    auto vals = split(line, '\t');
                    if (vals.size() != 2) {
                        std::cerr << "[odgi::heaps] line does not have a path.name and path.group value:"
                                  << std::endl << line << std::endl;
                        return 1;
                    }
                    auto& path_name = vals.front();
                    auto& group = vals.back();
                    if (!graph.has_path(vals.front())) {
                        std::cerr << "[odgi::heaps] no path '" << path_name << "' in graph" << std::endl;
                        return 1;
                    }
                    path_groups_map[group].push_back(graph.get_path_handle(path_name));
                }
            }
        } else if (group_by_haplotype) {
            graph.for_each_path_handle([&](const path_handle_t& p) {
                auto path_name = graph.get_path_name(p);
                // split and decide
                auto vals = split(path_name, '#');
                if (vals.size() == 1) { // Assume ctg
                    path_groups_map[vals.front()].push_back(p);
                } else if (vals.size() == 2) { // Assume sample#ctg
                    path_groups_map[vals.front()].push_back(p);
                } else if (vals.size() == 3) { // Assume sample#hap#ctg
                    auto key = vals[0] + vals[1];
                    path_groups_map[key].push_back(p);
                }
            });
        } else if (group_by_sample) {
            graph.for_each_path_handle([&](const path_handle_t& p) {
                auto path_name = graph.get_path_name(p);
                // split and decide
                auto vals = split(path_name, '#');
                path_groups_map[vals.front()].push_back(p);
            });
        } else {
            // no groups
            graph.for_each_path_handle([&](const path_handle_t& p) {
                path_groups_map[graph.get_path_name(p)].push_back(p);
            });
        }
        path_groups.reserve(path_groups_map.size());
        for (auto& g : path_groups_map) {
            path_groups.push_back(g.second);
        }
    }

    ska::flat_hash_map<path_handle_t, std::vector<interval_t>> intervals;
    if (_bed_targets) {
        std::ifstream bed(args::get(_bed_targets).c_str());
        std::string line;
        while (std::getline(bed, line)) {
            if (!line.empty()) {
                auto vals = split(line, '\t');
                if (vals.size() < 3) {
                    std::cerr << "[odgi::heaps]"
                              << "BED line does not have enough fields to define an interval"
                              << std::endl << line << std::endl;
                    return 1;
                }
                auto& path_name = vals[0];
                uint64_t start = std::stoul(vals[1]);
                uint64_t end = std::stoul(vals[2]);
                if (!graph.has_path(path_name)) {
                    //std::cerr << "[odgi::heaps] warning: no path '" << path_name << "' in graph" << std::endl;
                } else {
                    auto path = graph.get_path_handle(path_name);
                    intervals[path].push_back(interval_t(start, end));
                }
            }
        }
        // ensure sorted input
        std::vector<std::vector<interval_t>*> v;
        for (auto& i : intervals) {
            v.push_back(&i.second);
        }
#pragma omp parallel for
        for (auto& ivals : v) {
            std::sort(ivals->begin(), ivals->end());
        }
    }

    graph.set_number_of_threads(num_threads);

    std::cout << "permutation\tnth.genome\tbase.pairs" << std::endl;
    auto handle_output = [&](const std::vector<uint64_t>& vals, uint64_t perm_id) {
        int i = 0;
#pragma omp critical (cout)
        for (auto& v : vals) {
            std::cout << perm_id << "\t" << ++i << "\t" << v << std::endl;
        }
    };

    algorithms::for_each_heap_permutation(graph, path_groups, intervals, n_permutations, min_node_depth, handle_output);

    return 0;
}

static Subcommand odgi_heaps("heaps", "Path pangenome coverage permutations.",
                             PIPELINE, 3, main_heaps);


}
