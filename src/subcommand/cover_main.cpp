#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/cover.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_cover(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi cover";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("find a path cover of the graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the graph with the generated paths in this file", {'o', "out"});
    args::ValueFlag<uint64_t> num_paths_per_component(parser, "N", "number of paths to generate per component", {'n', "num-paths-per-component"});
    args::ValueFlag<uint64_t> node_window_size(parser, "N", "size of the node window to check each time a new path is extended (it has to be greater than or equal to 2)", {'k', "node-window-size"});
    args::Flag debug(parser, "debug", "Print information about the components and the progress to stdout", {'d', "debug"});

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
        std::cerr << "[odgi cover] error: Please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (!dg_out_file) {
        std::cerr << "[odgi cover] error: Please specify an output file to where to store the graph via -o=[FILE], --out=[FILE]." << std::endl;
        return 1;
    }

    uint64_t _num_paths_per_component = args::get(num_paths_per_component) ? args::get(num_paths_per_component) : algorithms::PATH_COVER_DEFAULT_N;
    uint64_t _node_window_size = args::get(node_window_size) ? args::get(node_window_size) : algorithms::PATH_COVER_DEFAULT_K;

    if (_node_window_size < 2){
        std::cerr << "[odgi cover] error: Please specify a node windows size greater than or equal to 2 -k=[N], --node-window-size=[N]." << std::endl;
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

    algorithms::path_cover(graph, _num_paths_per_component, _node_window_size, args::get(debug));
    
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

static Subcommand odgi_cover("cover", "find a path cover of the graph",
                            PIPELINE, 3, main_cover);


}
