#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "algorithms/prune.hpp"
#include "algorithms/coverage.hpp"
#include "algorithms/remove_high_degree.hpp"
#include "algorithms/cut_tips.hpp"

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
    args::Flag drop_paths(parser, "bool", "remove the paths from the graph", {'D', "drop-paths"});
    args::Flag cut_tips(parser, "bool", "remove nodes which are graph tips", {'T', "cut-tips"});
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
    assert(args::get(kmer_length));

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

    if (args::get(threads)) {
        omp_set_num_threads(args::get(threads));
    }
    if (args::get(max_degree)) {
        graph.clear_paths();
        algorithms::remove_high_degree_nodes(graph, args::get(max_degree));
    }
    if (args::get(max_furcations)) {
        std::vector<edge_t> to_prune = algorithms::find_edges_to_prune(graph, args::get(kmer_length), args::get(max_furcations));
        std::cerr << "edges to prune: " << to_prune.size() << std::endl;
        for (auto& edge : to_prune) {
            graph.destroy_edge(edge);
        }
        // we're just removing edges, so paths shouldn't be damaged
        std::cerr << "done prune" << std::endl;
    }
    if (args::get(min_coverage) || args::get(max_coverage)) {
        std::vector<handle_t> to_drop = algorithms::find_handles_exceeding_coverage_limits(graph, args::get(min_coverage), args::get(max_coverage));
        // remove the paths, because it's likely we have damaged some
        // and at present, we have no mechanism to reconstruct them
        graph.clear_paths();
        //std::cerr << "got " << to_drop.size() << " handles to drop" << std::endl;
        for (auto& handle : to_drop) {
            graph.destroy_handle(handle);
        }
    }
    if (args::get(cut_tips)) {
        algorithms::cut_tips(graph);
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
