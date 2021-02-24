#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/matrix_writer.hpp"

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

    if (!dg_in_file) {
        std::cerr << "[odgi::matrix] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
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

    algorithms::write_as_sparse_matrix(std::cout, graph, args::get(weight_by_edge_depth), args::get(weight_by_edge_delta));

    return 0;
}

static Subcommand odgi_matrix("matrix", "graph topology in sparse matrix form",
                              PIPELINE, 3, main_matrix);


}
