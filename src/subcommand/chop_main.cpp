#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "algorithms/simple_components.hpp"

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
    
    args::ArgumentParser parser("merge unbranching linear components into single nodes (drops paths)");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the graph self index in this file", {'o', "out"});
    args::ValueFlag<uint64_t> chop_to(parser, "N", "divide nodes to be shorter than this length", {'c', "chop-to"});
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

    graph_t graph;
    assert(argc > 0);
    assert(args::get(kmer_length));
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

    uint64_t max_node_length = args::get(chop_to);
    std::vector<handle_t> to_chop;
    graph.for_each_handle([&](const handle_t& handle) {
            if (graph.get_length(handle) > max_node_length) {
                to_chop.push_back(handle);
            }
        });

    for (auto& handle : to_chop) {
        // get divide points
        uint64_t length = graph.get_length(handle);
        std::vector<size_t> offsets;
        std::cerr << "node " << graph.get_id(handle) << " " << graph.get_length(handle) << std::endl;
        for (uint64_t i = max_node_length; i < length; i+=max_node_length) {
            std::cerr << "cutting at " << i << std::endl;
            offsets.push_back(i);
        }
        graph.divide_handle(handle, offsets);
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

static Subcommand odgi_chop("chop", "chop long nodes into short ones while preserving topology",
                            PIPELINE, 3, main_chop);


}
