#include "subcommand.hpp"
#include "graph.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "algorithms/simple_components.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_simplify(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi simplify";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("merge unbranching linear components into single nodes (drops paths)");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the graph self index in this file", {'o', "out"});
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

    graph_t graph;
    assert(argc > 0);
    assert(args::get(kmer_length));
    std::string infile = args::get(dg_in_file);
    if (infile.size()) {
        ifstream f(infile.c_str());
        graph.load(f);
        f.close();
    }
    if (args::get(threads)) {
        omp_set_num_threads(args::get(threads));
    }
    graph.clear_paths();
    std::vector<std::vector<handle_t>> linear_components = algorithms::simple_components(graph, 2);
    //std::cerr << "there are " << linear_components.size() << " components" << std::endl;
    uint64_t i = 0;
    for (auto& v : linear_components) {
        graph.combine_handles(v);
    }
    std::string outfile = args::get(dg_out_file);
    if (outfile.size()) {
        ofstream f(outfile.c_str());
        graph.serialize(f);
        f.close();
    }
    return 0;
}

static Subcommand odgi_simplify("simplify", "merge unbranching linear components into single nodes",
                                PIPELINE, 3, main_simplify);


}
