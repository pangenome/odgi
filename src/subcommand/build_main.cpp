#include "subcommand.hpp"
#include "odgi.hpp"
#include "gfa_to_handle.hpp"
#include "args.hxx"
#include <cstdio>
#include <algorithm>
#include <filesystem>
#include "algorithms/topological_sort.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_build(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi build";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("Construct a dynamic succinct variation graph in ODGI format from a GFAv1.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> gfa_file(mandatory_opts, "FILE", "GFAv1 FILE containing the nodes, edges and "
                                                          "paths to build a dynamic succinct variation graph from.", {'g', "gfa"});
    args::ValueFlag<std::string> dg_out_file(mandatory_opts, "FILE", "Write the dynamic succinct variation graph to this *FILE*. A file ending with *.og* is recommended.", {'o', "out"});
    args::Group graph_sorting(parser, "[ Graph Sorting ]");
    args::Flag optimize(graph_sorting, "optimize", "Compact the graph id space into a dense integer range.", {'O', "optimize"});
    args::Flag toposort(graph_sorting, "sort", "Apply a general topological sort to the graph and order the node ids"
                                        "  accordingly. A bidirected adaptation of Kahn’s topological sort (1962)"
                                        "  is used, which can handle components with no heads or tails. Here, both heads and tails are taken into account.", {'s', "sort"});
    args::Group threading(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> nthreads(threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
    args::Group processing_information(parser, "[ Processing Information ]");
    args::Flag progress(processing_information, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Flag debug(processing_information, "debug", "Verbosely print graph information to stderr. This includes the maximum"
                                                      "  node_id, the minimum node_id, the handle to node_id mapping, the"
                                                      "  deleted nodes and the path metadata.", {'d', "debug"});
    args::Group program_information(parser, "[ Program Information ]");
    args::HelpFlag help(program_information, "help", "Print a help message for odgi build.", {'h', "help"});
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
    
    //make_graph();
    assert(argc > 0);
    if (!gfa_file) {
		std::cerr << "[odgi::build] error: please specify an input file to load the graph from via -g=[FILE], --gfa=[FILE]." << std::endl;
		return 1;
    }
    if (!dg_out_file) {
        std::cerr << "[odgi::build] error: please specify an output file to store the graph via -o=[FILE], --out=[FILE]." << std::endl;
        return 1;
    }
    {
        const std::string gfa_filename = args::get(gfa_file);
        if (!std::filesystem::exists(gfa_filename)) {
            std::cerr << "[odgi::build] error: the given file \"" << gfa_filename << "\" does not exist. Please specify an existing input file via -g=[FILE], --gfa=[FILE]." << std::endl;
            return 1;
        }
        if (!gfa_filename.empty()) {
            gfa_to_handle(gfa_filename, &graph, args::get(optimize), args::get(nthreads), args::get(progress));
        }
    }

    const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;
    graph.set_number_of_threads(num_threads);

    if (args::get(toposort)) {
        graph.apply_ordering(algorithms::topological_order(&graph, true, args::get(progress)), true);
    }
    // here we should measure memory usage etc.
    if (args::get(debug)) {
        graph.display();
    }
    const std::string outfile = args::get(dg_out_file);
    if (!outfile.empty()) {
        if (outfile == "-") {
            graph.serialize(std::cout);
        } else {
            ofstream f(outfile.c_str());
            graph.serialize(f);
            f.close();
        }
    }
    return 0;
}

static Subcommand odgi_build("build", "Construct a dynamic succinct variation graph in ODGI format from a GFAv1.",
                              PIPELINE, 3, main_build);


}
