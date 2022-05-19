#include "subcommand.hpp"
#include <iostream>
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/break_cycles.hpp"
#include "utils.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_break(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    const std::string prog_name = "odgi break";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("Break cycles in the graph and drop its paths.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> odgi_in_file(mandatory_opts, "FILE", "Load the succinct variation graph "
                                                                      "in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::ValueFlag<std::string> odgi_out_file(mandatory_opts, "FILE", "Write the broken graph in ODGI format to FILE. A file ending of *.og* is recommended.", {'o', "out"});
    args::Group cycle_opts(parser, "[ Cycle Options ]");
    args::ValueFlag<uint64_t> max_cycle_size(cycle_opts, "N", "The maximum cycle length at which to break (default: 0).", {'c', "cycle-max-bp"});
    args::ValueFlag<uint64_t> max_search_bp(cycle_opts, "N", "The maximum search space of each BFS given in number of base pairs (default: 0).", {'s', "max-search-bp"});
    args::ValueFlag<uint64_t> repeat_up_to(cycle_opts, "N", "Iterate cycle breaking up to N times, or stop if no new edges are removed.", {'u', "repeat-up-to"});
    args::Flag show(cycle_opts, "show", "print edges we would remove", {'d', "show"});
	args::Group threading(parser, "[ Threading ]");
	args::ValueFlag<uint64_t> nthreads(threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
	args::Group processing_info_opts(parser, "[ Processing Information ]");
	args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi break.", {'h', "help"});
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

    if (!odgi_in_file) {
        std::cerr << "[odgi::break] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (!args::get(show) && !odgi_out_file) {
        std::cerr << "[odgi::break] error: please specify an output file to where to store the broken graph via -o=[FILE], --out=[FILE]." << std::endl;
        return 1;
    }

    if (args::get(show) && odgi_out_file) {
        std::cerr << "[odgi::break] error: please do not specify an output file when edges we would remove are printed (-d/--show)." << std::endl;
        return 1;
    }

	const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

    graph_t graph;
    assert(argc > 0);
    {
        const std::string infile = args::get(odgi_in_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
				utils::handle_gfa_odgi_input(infile, "break", args::get(progress), num_threads, graph);
            }
        }
    }

    const uint64_t iter_max = args::get(repeat_up_to) ? args::get(repeat_up_to) : 1;

    // break cycles
    if (args::get(show)) {
        std::vector<edge_t> cycle_edges
            = algorithms::edges_inducing_cycles(graph, args::get(max_cycle_size), args::get(max_search_bp));

        for (auto& e : cycle_edges) {
            std::cout << graph.get_id(e.first) << (graph.get_is_reverse(e.first)?"-":"+")
                      << " -> "
                      << graph.get_id(e.second) << (graph.get_is_reverse(e.second)?"-":"+")
                      << std::endl;
        }
    } else {
        const uint64_t removed_edges
            = algorithms::break_cycles(graph, args::get(max_cycle_size), args::get(max_search_bp), iter_max);
        if (removed_edges > 0) {
            graph.clear_paths();
        }

        {
            const std::string outfile = args::get(odgi_out_file);
            if (!outfile.empty()) {
                if (outfile == "-") {
                    graph.serialize(std::cout);
                } else {
                    ofstream f(outfile.c_str());
                    graph.serialize(f);
                    f.close();
                }
            }
        }
    }

    return 0;
}

static Subcommand odgi_break("break", "Break cycles in the graph and drop its paths.",
                              PIPELINE, 3, main_break);


}
