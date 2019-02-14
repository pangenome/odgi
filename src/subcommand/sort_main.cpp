#include "subcommand.hpp"
#include "graph.hpp"
#include "args.hxx"
#include "algorithms/topological_sort.hpp"

namespace dsgvg {

using namespace dsgvg::subcommand;

int main_sort(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "dsgvg sort";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("metrics describing variation graphs");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the graph in this file", {'o', "out"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::Flag show_sort(parser, "show", "write the sort order mapping", {'S', "show"});
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
        ifstream f(infile.c_str());
        graph.load(f);
        f.close();
    }
    if (args::get(show_sort)) {
        vector<handle_t> order = algorithms::topological_order(&graph);
        for (auto& handle : order) {
            std::cout << graph.get_id(handle) << std::endl;
        }
    }
    std::string outfile = args::get(dg_out_file);
    if (outfile.size()) {
        graph.apply_ordering(algorithms::topological_order(&graph), true);
        ofstream f(outfile.c_str());
        graph.serialize(f);
        f.close();
    }
    return 0;
}

static Subcommand dsgvg_build("sort", "topologically order the graph",
                              PIPELINE, 3, main_sort);


}
