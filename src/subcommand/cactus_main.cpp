#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "cactus.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_cactus(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi cactus";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("build the cactus graph transformation of the graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the index from this file", {'i', "idx"});
    //args::Flag to_gfa(parser, "to_gfa", "write the graph to stdout in GFA format", {'g', "to-gfa"});
    //args::Flag display(parser, "display", "show internal structures", {'d', "display"});

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
    auto g = handle_graph_to_cactus(graph, {});
    stCactusGraph* cactus_graph = g.first;
    stList* telomeres = g.second;
    // 
    stList_destruct(telomeres);
    stCactusGraph_destruct(cactus_graph);


    return 0;
}

static Subcommand odgi_cactus("cactus", "cactus transformation of the graph",
                              PIPELINE, 3, main_cactus);


}
