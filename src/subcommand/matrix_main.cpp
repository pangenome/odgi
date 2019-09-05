#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/bin_path_info.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_matrix(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi matrix";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("write the graph topology in sparse matrix formats");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the graph in this file", {'o', "out"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::Flag weight_by_edge_depth(parser, "edge-depth-weight", "weight edges by their path depth", {'e', "edge-depth-weight"});
    args::Flag weight_by_edge_delta(parser, "delta-weight", "weight edges by the inverse id delta", {'d', "delta-weight"});
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

    uint64_t edge_count = 0;
    graph.for_each_edge([&](const edge_t& edge) {
            ++edge_count;
        });
    std::cout << graph.max_node_id() << " " << graph.max_node_id() << " " << edge_count*2 << std::endl;
    graph.for_each_edge([&](const edge_t& edge) {
            // how many paths cross the edge?
            double weight = 0;
            if (args::get(weight_by_edge_depth)) {
                graph.for_each_step_on_handle(edge.first, [&](const step_handle_t& step) {
                        if (graph.get_handle_of_step(step) == edge.first
                            && graph.has_next_step(step)) {
                            handle_t next = graph.get_handle_of_step(graph.get_next_step(step));
                            if (next == edge.second) {
                                ++weight;
                            }
                        } else if (graph.get_handle_of_step(step) == graph.flip(edge.first)
                                   && graph.has_previous_step(step)) {
                            handle_t prev = graph.get_handle_of_step(graph.get_previous_step(step));
                            if (prev == graph.flip(edge.second)) {
                                ++weight;
                            }
                        }
                    });
            } else {
                weight = 1;
            }
            if (args::get(weight_by_edge_delta)) {
                double delta = std::abs(graph.get_id(edge.first) - graph.get_id(edge.second));
                if (delta == 0) delta = 1;
                weight = 1 / delta;
            }
            std::cout << graph.get_id(edge.first) << " " << graph.get_id(edge.second) << " " << weight << std::endl;
            std::cout << graph.get_id(edge.second) << " " << graph.get_id(edge.first) << " " << weight << std::endl;
        });

    return 0;
}

static Subcommand odgi_matrix("matrix", "graph topology in sparse matrix form",
                              PIPELINE, 3, main_matrix);


}
