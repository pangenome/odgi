#include "subcommand.hpp"
#include "graph.hpp"
#include "algorithms/kmer.hpp"
#include "args.hxx"

namespace dsgvg {

using namespace dsgvg::subcommand;

int main_kmers(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "dsgvg kmers";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("show and characterize the kmer space of the graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<uint64_t> kmer_length(parser, "kmer_length", "the length of the kmers to generate", {'k', "kmer-length"});
    args::Flag kmers_stdout(parser, "kmers_stdout", "write the kmers to stdout", {'c', "stdout"});

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
    assert(args::get(kmer_length));
    std::string infile = args::get(dg_in_file);
    if (infile.size()) {
        ifstream f(infile.c_str());
        graph.load(f);
        f.close();
    }
    if (args::get(kmers_stdout)) {
        algorithms::for_each_kmer(graph, args::get(kmer_length), [&](const kmer_t& kmer) {
#pragma omp critical
                std::cout << kmer << std::endl;
            });
    }
    return 0;
}

static Subcommand dsgvg_kmers("kmers", "process and dump the kmers of the graph",
                              PIPELINE, 3, main_kmers);


}
