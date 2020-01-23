#include "subcommand.hpp"
#include <iostream>
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/break_cycles.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_break(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi break";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("break cycles in the graph (drops paths)");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> odgi_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> odgi_out_file(parser, "FILE", "store the graph in this file", {'o', "out"});
    args::ValueFlag<uint64_t> max_cycle_size(parser, "N", "maximum cycle length to break", {'c', "cycle-max-bp"});
    args::ValueFlag<uint64_t> max_search_bp(parser, "N", "maximum number of bp per BFS from any node", {'s', "max-search-bp"});
    args::ValueFlag<uint64_t> repeat_up_to(parser, "N", "iterate cycle breaking up to N times, or stop if no new edges are removed", {'u', "repeat-up-to"});
    args::Flag show(parser, "show", "print edges we would remove", {'d', "show"});

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
    std::string infile = args::get(odgi_in_file);
    if (infile.size()) {
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
            ifstream f(infile.c_str());
            graph.deserialize(f);
            f.close();
        }
    }

    uint64_t iter_max = args::get(repeat_up_to) ? args::get(repeat_up_to) : 1;

    // break cycles
    if (args::get(show)) {
        std::vector<edge_t> cycle_edges
            = algorithms::edges_inducing_cycles(graph,
                                                args::get(max_cycle_size),
                                                args::get(max_search_bp));
        for (auto& e : cycle_edges) {
            std::cout << graph.get_id(e.first) << (graph.get_is_reverse(e.first)?"-":"+")
                      << " -> "
                      << graph.get_id(e.second) << (graph.get_is_reverse(e.second)?"-":"+")
                      << std::endl;
        }
    } else {
        uint64_t removed_edges
            = algorithms::break_cycles(graph,
                                       args::get(max_cycle_size),
                                       args::get(max_search_bp),
                                       iter_max);
        if (removed_edges > 0) {
            graph.clear_paths();
        }
    }

    std::string outfile = args::get(odgi_out_file);
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

static Subcommand odgi_break("break", "break cycles in the graph",
                              PIPELINE, 3, main_break);


}
