#include "subcommand.hpp"
#include "graph.hpp"
//#include "gfakluge.hpp"
#include "args.hxx"
//#include "io_helper.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_view(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi view";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("projection of graphs into other formats");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    //args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the index in this file", {'o', "out"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the index from this file", {'i', "idx"});
    //args::ValueFlag<std::string> seqs(parser, "FILE", "the sequences used to generate the alignments", {'s', "seqs"});
    //args::ValueFlag<std::string> base(parser, "FILE", "build graph using this basename", {'b', "base"});
    //args::ValueFlag<uint64_t> num_threads(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    //args::ValueFlag<uint64_t> repeat_max(parser, "N", "limit transitive closure to include no more than N copies of a given input base", {'r', "repeat-max"});
    //args::ValueFlag<uint64_t> aln_keep_n_longest(parser, "N", "keep up to the N-longest alignments overlapping each query position", {'k', "aln-keep-n-longest"});
    //args::ValueFlag<uint64_t> aln_min_length(parser, "N", "ignore alignments shorter than this", {'m', "aln-min-length"});
    args::Flag to_gfa(parser, "to_gfa", "write the graph to stdout in GFA format", {'g', "to-gfa"});
    //args::Flag summarize(parser, "summarize", "summarize the graph properties and dimensions", {'S', "summarize"});
    args::Flag display(parser, "display", "show internal structures", {'d', "display"});
    //args::Flag progress(parser, "progress", "show progress updates", {'p', "progress"});
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
    if (args::get(display)) {
        graph.display();
    }
    if (args::get(to_gfa)) {
        graph.to_gfa(std::cout);
    }

    return 0;
}

static Subcommand odgi_view("view", "projection of graphs into other formats",
                             PIPELINE, 3, main_view);


}
