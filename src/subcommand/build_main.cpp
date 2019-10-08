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
    
    args::ArgumentParser parser("construct a dynamic succinct variation graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> gfa_file(parser, "FILE", "construct the graph from this GFA input file", {'g', "gfa"});
    args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the graph self index in this file", {'o', "out"});
    args::Flag to_gfa(parser, "to_gfa", "write the graph to stdout in GFA format", {'G', "to-gfa"});
    args::Flag toposort(parser, "sort", "apply generalized topological sort to the graph and set node ids to order", {'s', "sort"});
    args::Flag debug(parser, "debug", "enable debugging", {'d', "debug"});
    args::Flag progress(parser, "progress", "show progress updates", {'p', "progress"});
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

    /*
    size_t n_threads = args::get(num_threads);
    if (n_threads) {
        omp_set_num_threads(args::get(num_threads));
    } else {
        omp_set_num_threads(1);
    }
    */
    graph_t graph;
    
    //make_graph();
    assert(argc > 0);
    std::string gfa_filename = args::get(gfa_file);
    if (gfa_filename.size()) {
        gfa_to_handle(gfa_filename, &graph, args::get(progress));
    }
    if (args::get(progress)) {
        std::cerr << std::endl;
    }
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
