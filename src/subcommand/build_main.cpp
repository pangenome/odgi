#include "subcommand.hpp"
#include "odgi.hpp"
#include "gfa_to_handle.hpp"
#include "args.hxx"
#include <cstdio>
#include <algorithm>
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
    args::Group mandatory(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> gfa_file(mandatory, "FILE", "GFAv1 FILE containing the nodes, edges and "
                                                          "paths to build a dynamic succinct variation graph from.", {'g', "gfa"});
    args::ValueFlag<std::string> dg_out_file(mandatory, "FILE", "Write the dynamic succinct variation graph to this file. A file ending with *.og* is recommended.", {'o', "out"});
    args::Group graph_files_io(parser, "[ Graph Files IO ]");
    args::Flag to_gfa(graph_files_io, "to_gfa", "Write the graph to stdout in GFAv1 format.", {'G', "to-gfa"});
    args::Group graph_sorting(parser, "[ Graph Sorting ]");
    args::Flag toposort(graph_sorting, "sort", "Apply a general topological sort to the graph and order the node ids"
                                        "  accordingly. A bidirected adaptation of Kahn’s topological sort (1962)"
                                        "  is used, which can handle components with no heads or tails. Here, both heads and tails are taken into account.", {'s', "sort"});
    args::Group threading(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> nthreads(threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
    args::Group program_information(parser, "[ Program Information ]");
    args::Flag progress(program_information, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::HelpFlag help(program_information, "help", "Display this help summary.", {'h', "help"});
    args::Group processing_information(parser, "[ Processing Information ]");
    args::Flag debug(processing_information, "debug", "Verbosely print graph information to stderr. This includes the maximum"
                                                      "  node_id, the minimum node_id, the handle to node_id mapping, the"
                                                      "  deleted nodes and the path metadata.", {'d', "debug"});
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
    if (!dg_out_file) {
        std::cerr << "Please specify an output file to store the graph via -o=[FILE], --out=[FILE]." << std::endl;
        return 1;
    }
    std::string gfa_filename = args::get(gfa_file);
    if (gfa_filename.size()) {
        gfa_to_handle(gfa_filename, &graph, args::get(nthreads), args::get(progress));
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
    if (args::get(to_gfa)) {
        graph.to_gfa(std::cout);
    }
    std::string outfile = args::get(dg_out_file);
    if (outfile.size()) {
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

static Subcommand odgi_build("build", "build dynamic succinct variation graph",
                              PIPELINE, 3, main_build);


}
