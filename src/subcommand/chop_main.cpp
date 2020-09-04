#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "algorithms/chop.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_chop(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi chop";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("divide nodes into smaller pieces");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the graph self index in this file", {'o', "out"});
    args::ValueFlag<uint64_t> chop_to(parser, "N", "divide nodes to be shorter than this length", {'c', "chop-to"});
    args::ValueFlag<uint64_t> nthreads(parser, "N", "number of threads to use for the parallel sorter", {'t', "threads"});
    args::Flag debug(parser, "debug", "print information about the components", {'d', "debug"});

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
        std::cerr << "[odgi chop] error: Please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (!dg_out_file) {
        std::cerr << "[odgi chop] error: Please specify an output file to where to store the graph via -o=[FILE], --out=[FILE]." << std::endl;
        return 1;
    }

    if (!chop_to) {
        std::cerr << "[odgi chop] error: Please specify a node chop length via -c=[N], --chop-to=[N]." << std::endl;
        return 1;
    }

    graph_t graph;
    assert(argc > 0);
    std::string infile = args::get(dg_in_file);
    if (infile.size()) {
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
            ifstream f(infile.c_str());
            graph.deserialize(f);
            f.close();
        }
    }

    uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

    algorithms::chop(graph, args::get(chop_to), num_threads, args::get(debug));
    
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

static Subcommand odgi_chop("chop", "chop long nodes into short ones while preserving topology and node order",
                            PIPELINE, 3, main_chop);


}
