#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "split.hpp"
#include <omp.h>

#include "src/algorithms/subgraph/extract.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_degree(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi degree";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("Describe the graph in terms of node degree.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> og_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "input"});
    args::Group summary_opts(parser, "[ Summary Options ]");
    args::Flag _summarize(summary_opts, "summarize", "Summarize the graph properties and dimensions. Print to stdout the node.id and the node.degree.", {'S', "summarize"});
    args::ValueFlag<std::string> _windows_in(summary_opts, "LEN:MIN:MAX",
                                             "Print to stdout a BED file of path intervals where the degree is between MIN and MAX, "
                                             "merging regions not separated by more than LEN bp.",
                                             {'w', "windows-in"});
    args::ValueFlag<std::string> _windows_out(summary_opts, "LEN:MIN:MAX",
                                              "Print to stdout a BED file of path intervals where the degree is outside of MIN and MAX, "
                                              "merging regions not separated by more than LEN bp.",
                                              {'W', "windows-out"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> _num_threads(threading_opts, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
	args::Group processing_info_opts(parser, "[ Processing Information ]");
	args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi degree.", {'h', "help"});
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

    if (!og_file) {
        std::cerr << "[odgi::degree] error: please specify a target graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (_windows_in && _windows_out) {
        std::cerr << "[odgi::degree] error: please specify -w/--windows-in or -W/--windows-out, not both." << std::endl;
        return 1;
    }

    if (_summarize && (_windows_in || _windows_out)) {
        std::cerr << "[odgi::degree] error: please specify -S/--summarize without specifying windows-in or -W/--windows-out." << std::endl;
        return 1;
    }

    uint64_t windows_in_len = 0, windows_in_min = 0, windows_in_max = 0;
    if (_windows_in) {
        if (!algorithms::check_and_get_windows_in_out_parameter(args::get(_windows_in), windows_in_len, windows_in_min, windows_in_max)) {
            std::cerr << "[odgi::degree] error: please specify a valid string (LEN:MIN:MAX) for the -w/--windows-in option." << std::endl;
            return 1;
        }
    }

    uint64_t windows_out_len = 0, windows_out_min = 0, windows_out_max = 0;
    if (_windows_out) {
        if (!algorithms::check_and_get_windows_in_out_parameter(args::get(_windows_out), windows_out_len, windows_out_min, windows_out_max)) {
            std::cerr << "[odgi::degree] error: please specify a valid string (LEN:MIN:MAX) for the -W/--windows-out option." << std::endl;
            return 1;
        }
    }

	const uint64_t num_threads = args::get(_num_threads) ? args::get(_num_threads) : 1;

	odgi::graph_t graph;
    assert(argc > 0);
    if (!args::get(og_file).empty()) {
        std::string infile = args::get(og_file);
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
			utils::handle_gfa_odgi_input(infile, "degree", args::get(progress), num_threads, graph);
        }
    }

    omp_set_num_threads((int) num_threads);

    if (_summarize) {
        uint64_t total_edges = 0;
        uint64_t min_degree = std::numeric_limits<uint64_t>::max();
        uint64_t max_degree = std::numeric_limits<uint64_t>::min();
        graph.for_each_handle(
                [&](const handle_t& handle) {
                    uint64_t degree = graph.get_degree(handle, false) + graph.get_degree(handle, true);
                    total_edges += degree;
                    min_degree = std::min(min_degree, degree);
                    max_degree = std::max(max_degree, degree);
                });
        std::cout << "#node.count\tedge.count\tavg.degree\tmin.degree\tmax.degree" << std::endl
                  << graph.get_node_count() << "\t"
                  << total_edges / 2 << "\t" // we double-count edges
                  << (double) total_edges / (double)graph.get_node_count() << "\t"
                  << min_degree << "\t"
                  << max_degree
                  << std::endl;
    } else  if (_windows_in || _windows_out) {
        std::vector<path_handle_t> paths;
        paths.reserve(graph.get_path_count());
        graph.for_each_path_handle([&](const path_handle_t path) {
            paths.push_back(path);
        });

        // precompute degrees for all handles in parallel
        std::vector<uint64_t> degrees(graph.get_node_count() + 1);
        graph.for_each_handle(
            [&](const handle_t& h) {
                auto id = graph.get_id(h);
                if (id >= degrees.size()) {
                	// require optimized graph to use vector rather than a hash table
                	std::cerr << "[odgi::degree] error: graph is not optimized, apply 'odgi sort' with -O, --optimize." << std::endl;
                	assert(false);
                }
                degrees[id] = graph.get_degree(h, false) + graph.get_degree(h, true);
            }, true);

        auto in_bounds =
            [&](const handle_t &handle) {
                uint64_t degree = degrees[graph.get_id(handle)];
                return _windows_in ? (degree >= windows_in_min && degree <= windows_in_max) : (degree < windows_out_min || degree > windows_out_max);
            };

        std::cout << "#path\tstart\tend" << std::endl;

        algorithms::windows_in_out(graph, paths, in_bounds, _windows_in ? windows_in_len : windows_out_len,
                                   [&](const std::vector<path_range_t>& path_ranges) {
#pragma omp critical (cout)
                                       for (auto path_range : path_ranges) {
                                           std::cout << graph.get_path_name(path_range.begin.path) << "\t"
                                                     << path_range.begin.offset << "\t"
                                                     << path_range.end.offset << std::endl;
                                       }
                                   }, num_threads);
    } else {
        std::cout << "#node.id\tnode.degree" << std::endl;
        graph.for_each_handle(
                [&](const handle_t& handle) {
#pragma omp critical (cout)
                    std::cout << graph.get_id(handle) << "\t"
                              << graph.get_degree(handle, false) + graph.get_degree(handle, true)
                              << std::endl;
                }, num_threads > 1);
    }

    return 0;
}

static Subcommand odgi_degree("degree", "Describe the graph in terms of node degree.",
                              PIPELINE, 3, main_degree);


}
