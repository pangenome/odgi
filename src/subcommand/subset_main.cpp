#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>

namespace odgi {

using namespace odgi::subcommand;

int main_subset(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi subset";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("extract subsets of the graph as defined by query criteria");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the graph self index in this file", {'o', "out"});
    args::ValueFlag<std::string> node_list(parser, "FILE", "a file with one node id per line", {'l', "node-list"});
    args::ValueFlag<uint64_t> target_node(parser, "ID", "a single node from which to begin our traversal", {'n', "node"});
    args::ValueFlag<uint64_t> context_size(parser, "N", "the number of steps away from our initial subgraph that we should collect", {'c', "context"});
    args::ValueFlag<uint64_t> threads(parser, "N", "number of threads to use", {'t', "threads"});

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
        std::cerr << "[odgi subset] error: Please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (!dg_out_file) {
        std::cerr << "[odgi subset] error: Please specify an output file to where to store the unchopped graph via -o=[FILE], --out=[FILE]." << std::endl;
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
    if (args::get(threads)) {
        omp_set_num_threads(args::get(threads));
    }

    graph_t subset;
    // collect the new graph
    if (args::get(target_node)) {
        uint64_t node_id = args::get(target_node);
        handle_t handle = graph.get_handle(node_id);
        subset.create_handle(graph.get_sequence(handle), node_id);
    }
    if (!args::get(node_list).empty()) {
        ifstream nodes(args::get(node_list));
        std::string next;
        while (std::getline(nodes, next)) {
            uint64_t id = std::stol(next);
            if (graph.has_node(id)) {
                handle_t h = graph.get_handle(id);
                subset.create_handle(graph.get_sequence(h), id);
            } else {
                std::cerr << "[odgi subset] Warning, cannot find node " << id << std::endl;
            }
        }
    }
    if (args::get(context_size)) {
        uint64_t context = args::get(context_size);
        for (uint64_t i = 0; i < context; ++i) {
            // get the edges and connected nodes from the graph to fill out th subgraph
            // for each subset node
            // get its component in the old graph
            // and add it to the subset
            std::vector<handle_t> curr_handles;
            subset.for_each_handle([&](const handle_t& h) {
                    curr_handles.push_back(h);
                });
            for (auto& h : curr_handles) {
                handle_t old_h = graph.get_handle(subset.get_id(h));
                graph.follow_edges(old_h, false, [&](const handle_t& c) {
                        handle_t x = subset.create_handle(graph.get_sequence(c), graph.get_id(c));
                        subset.create_edge(h, x);
                    });
                graph.follow_edges(old_h, true, [&](const handle_t& c) {
                        handle_t x = subset.create_handle(graph.get_sequence(c), graph.get_id(c));
                        subset.create_edge(x, h);
                    });
            }
        }
    }

    std::string outfile = args::get(dg_out_file);
    if (outfile.size()) {
        if (outfile == "-") {
            subset.serialize(std::cout);
        } else {
            ofstream f(outfile.c_str());
            subset.serialize(f);
            f.close();
        }
    }
    return 0;
}

static Subcommand odgi_subset("subset", "extract subsets of the graph as defined by query criteria",
                              PIPELINE, 3, main_subset);


}
