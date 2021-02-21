#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/prune.hpp"
#include "algorithms/coverage.hpp"
#include "algorithms/remove_high_degree.hpp"
#include "algorithms/cut_tips.hpp"
#include "algorithms/remove_isolated.hpp"
#include "algorithms/expand_context.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_prune(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi prune";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("remove complex parts of the graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the graph self index in this file", {'o', "out"});
    args::ValueFlag<uint64_t> kmer_length(parser, "K", "the length of the kmers to consider", {'k', "kmer-length"});
    args::ValueFlag<uint64_t> max_furcations(parser, "N", "break at edges that would be induce this many furcations in a kmer", {'e', "max-furcations"});
    args::ValueFlag<uint64_t> max_degree(parser, "N", "remove nodes that have degree greater that this level", {'d', "max-degree"});
    args::ValueFlag<uint64_t> min_coverage(parser, "N", "remove nodes covered by fewer than this number of path steps", {'c', "min-coverage"});
    args::ValueFlag<uint64_t> max_coverage(parser, "N", "remove nodes covered by more than this number of path steps", {'C', "max-coverage"});
    args::Flag edge_coverage(parser, "bool", "remove edges outside of the min and max coverage (rather than nodes)", {'E', "edge-coverage"});
    args::ValueFlag<uint64_t> best_edges(parser, "N", "keep only the N most-covered inbound and outbound edge of each node", {'b', "best-edges"});
    args::ValueFlag<uint64_t> expand_steps(parser, "N", "include nodes within this many steps of a component passing the prune thresholds", {'s', "expand-steps"});
    args::ValueFlag<uint64_t> expand_length(parser, "N", "include nodes within this graph distance of a component passing the prune thresholds", {'l', "expand-length"});
    args::ValueFlag<uint64_t> expand_path_length(parser, "N", "include nodes within this path length of a component passing the prune thresholds", {'p', "expand-path-length"});
    args::Flag drop_paths(parser, "bool", "remove the paths from the graph", {'D', "drop-paths"});
    args::Flag cut_tips(parser, "bool", "remove nodes which are graph tips", {'T', "cut-tips"});
    args::ValueFlag<uint64_t> cut_tips_min_coverage(parser, "bool", "remove nodes which are graph tips and have less than this path coverage", {'m', "cut-tips-min-coverage"});
    args::Flag remove_isolated(parser, "bool", "remove isolated nodes covered by a single path", {'I', "remove-isolated"});
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

    assert(argc > 0);

    if (!dg_in_file) {
        std::cerr << "[odgi::prune] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (!dg_out_file) {
        std::cerr
                << "[odgi::prune] error: please specify an output file to where to store the pruned graph via -o=[FILE], --out=[FILE]."
                << std::endl;
        return 1;
    }

    graph_t graph;
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

    int n_threads = threads ? args::get(threads) : 1;
    omp_set_num_threads(n_threads);

    if (args::get(max_degree)) {
        graph.clear_paths();
        algorithms::remove_high_degree_nodes(graph, args::get(max_degree));
    }
    if (args::get(max_furcations)) {
        std::vector<edge_t> to_prune = algorithms::find_edges_to_prune(graph, args::get(kmer_length), args::get(max_furcations), n_threads);
        //std::cerr << "edges to prune: " << to_prune.size() << std::endl;
        for (auto& edge : to_prune) {
            graph.destroy_edge(edge);
        }
        // we're just removing edges, so paths shouldn't be damaged
        //std::cerr << "done prune" << std::endl;
    }
    if (args::get(min_coverage) || args::get(max_coverage) || args::get(best_edges)) {
        std::vector<handle_t> handles_to_drop;
        std::vector<edge_t> edges_to_drop_coverage, edges_to_drop_best;
        if (args::get(edge_coverage)) {
            edges_to_drop_coverage = algorithms::find_edges_exceeding_coverage_limits(graph, args::get(min_coverage), args::get(max_coverage));
        } else {
            handles_to_drop = algorithms::find_handles_exceeding_coverage_limits(graph, args::get(min_coverage), args::get(max_coverage));
        }
        if (args::get(best_edges)) {
            edges_to_drop_best = algorithms::keep_mutual_best_edges(graph, args::get(best_edges));
        }
        // TODO this needs fixing
        // we should split up the paths rather than drop them
        // remove the paths, because it's likely we have damaged some
        // and at present, we have no mechanism to reconstruct them
        auto do_destroy =
            [&](void) {
                if (args::get(min_coverage) == 1 && args::get(max_coverage) == 0) {
                    // we could not have damaged any paths
                } else {
                    graph.clear_paths();
                }
                //std::cerr << "got " << to_drop.size() << " handles to drop" << std::endl;
                for (auto& edge : edges_to_drop_coverage) {
                    graph.destroy_edge(edge);
                }
                for (auto& edge : edges_to_drop_best) {
                    graph.destroy_edge(edge);
                }
                for (auto& handle : handles_to_drop) {
                    graph.destroy_handle(handle);
                }
            };
        if (args::get(expand_steps)) {
            graph_t source;
            source.copy(graph);
            do_destroy();
            algorithms::expand_context(&source, &graph, args::get(expand_steps), true);
        } else if (args::get(expand_length)) {
            graph_t source;
            source.copy(graph);
            do_destroy();
            algorithms::expand_context(&source, &graph, args::get(expand_length), false);
        } else if (args::get(expand_path_length)) {
            graph_t source;
            source.copy(graph);
            do_destroy();
            algorithms::expand_context_with_paths(&source, &graph, args::get(expand_length), false);
        } else {
            do_destroy();
        }
    }
    if (args::get(cut_tips)) {
        algorithms::cut_tips(graph, args::get(cut_tips_min_coverage));
        graph.optimize();
    }
    if (args::get(remove_isolated)) {
        algorithms::remove_isolated_paths(graph);
        graph.optimize();
    }
    if (args::get(drop_paths)) {
        graph.clear_paths();
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

static Subcommand odgi_prune("prune", "prune the graph based on coverage or topological complexity",
                              PIPELINE, 3, main_prune);


}
