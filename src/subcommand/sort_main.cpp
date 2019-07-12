#include "subcommand.hpp"
#include "odgi.hpp"
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
    
    args::ArgumentParser parser("variation graph sorts");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the graph in this file", {'o', "out"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::Flag show_sort(parser, "show", "write the sort order mapping", {'S', "show"});
    args::Flag cycle_breaking(parser, "cycle_breaking", "use a cycle breaking sort", {'b', "cycle-breaking"});
    args::Flag eades(parser, "eades", "use eades algorithm", {'e', "eades"});
    args::Flag lazy(parser, "lazy", "use lazy topological algorithm (DAG only)", {'l', "lazy"});
    args::Flag two(parser, "two", "use two-way (max of head-first and tail-first) topological algorithm", {'w', "two-way"});
    args::Flag paths_by_min_node_id(parser, "paths-min", "sort paths by their lowest contained node id", {'P', "paths-min"});
    args::Flag paths_by_max_node_id(parser, "paths-max", "sort paths by their highest contained node id", {'M', "paths-max"});
    args::Flag paths_by_avg_node_id(parser, "paths-avg", "sort paths by their average contained node id", {'A', "paths-avg"});
    args::Flag paths_by_avg_node_id_rev(parser, "paths-avg-rev", "sort paths in reverse by their average contained node id", {'R', "paths-avg-rev"});
    args::ValueFlag<std::string> path_delim(parser, "path-delim", "sort paths in bins by their prefix up to this delemiter", {'D', "path-delim"});
    args::Flag progress(parser, "progress", "display progress of the sort", {'p', "progress"});
    args::Flag optimize(parser, "optimize", "use the MutableHandleGraph::optimize method", {'O', "optimize"});
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
    if (args::get(show_sort)) {
        vector<handle_t> order = (args::get(lazy) ? algorithms::lazy_topological_order(&graph) : algorithms::topological_order(&graph));
        for (auto& handle : order) {
            std::cout << graph.get_id(handle) << std::endl;
        }
    }
    std::string outfile = args::get(dg_out_file);
    if (outfile.size()) {
        if (args::get(eades)) {
            graph.apply_ordering(algorithms::eades_algorithm(&graph), true);
        } else if (args::get(lazy)) {
            graph.apply_ordering(algorithms::lazy_topological_order(&graph), true);
        } else if (args::get(two)) {
            graph.apply_ordering(algorithms::two_way_topological_order(&graph), true);
        } else if (args::get(optimize)) {
            graph.optimize();
        } else {
            graph.apply_ordering(algorithms::topological_order(&graph, true, args::get(progress)), true);
        }
        if (args::get(cycle_breaking)) {
            graph.apply_ordering(algorithms::cycle_breaking_sort(graph), true);
        }
        if (args::get(paths_by_min_node_id)) {
            graph.apply_path_ordering(algorithms::prefix_and_id_ordered_paths(graph, args::get(path_delim), false, false));
        }
        if (args::get(paths_by_max_node_id)) {
            graph.apply_path_ordering(algorithms::prefix_and_id_ordered_paths(graph, args::get(path_delim), false, true));
        }
        if (args::get(paths_by_avg_node_id)) {
            graph.apply_path_ordering(algorithms::prefix_and_id_ordered_paths(graph, args::get(path_delim), true, false));
        }
        if (args::get(paths_by_avg_node_id_rev)) {
            graph.apply_path_ordering(algorithms::prefix_and_id_ordered_paths(graph, args::get(path_delim), true, true));
        }
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

static Subcommand odgi_build("sort", "topologically order the graph",
                              PIPELINE, 3, main_sort);


}
