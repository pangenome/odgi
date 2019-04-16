#include "subcommand.hpp"
#include "graph.hpp"
#include "args.hxx"
#include "algorithms/topological_sort.hpp"
#include "algorithms/eades_algorithm.hpp"
#include "algorithms/cycle_breaking_sort.hpp"
#include "algorithms/id_ordered_paths.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_sort(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi sort";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("metrics describing variation graphs");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the graph in this file", {'o', "out"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::Flag show_sort(parser, "show", "write the sort order mapping", {'S', "show"});
    args::Flag cycle_breaking(parser, "cycle_breaking", "use a cycle breaking sort", {'b', "cycle-breaking"});
    args::Flag eades(parser, "eades", "use eades algorithm", {'e', "eades"});
    args::Flag paths_by_min_node_id(parser, "paths-min", "sort paths by their lowest contained node id", {'P', "paths-min"});
    args::Flag paths_by_max_node_id(parser, "paths-max", "sort paths by their highest contained node id", {'M', "paths-max"});
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
        if (args::get(eades)) {
            graph.apply_ordering(algorithms::eades_algorithm(&graph), true);
        } else {
            graph.apply_ordering(algorithms::topological_order(&graph), true);
        }
        if (args::get(cycle_breaking)) {
            graph.apply_ordering(algorithms::cycle_breaking_sort(graph), true);
        }
        if (args::get(paths_by_min_node_id)) {
            graph.apply_path_ordering(algorithms::id_ordered_paths(graph, false));
        }
        if (args::get(paths_by_max_node_id)) {
            graph.apply_path_ordering(algorithms::id_ordered_paths(graph, true));
        }
        ofstream f(outfile.c_str());
        graph.serialize(f);
        f.close();
    }
    return 0;
}

static Subcommand odgi_build("sort", "topologically order the graph",
                              PIPELINE, 3, main_sort);


}
