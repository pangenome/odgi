#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "split.hpp"
#include <omp.h>

namespace odgi {

using namespace odgi::subcommand;

int main_degree(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi degree";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("describe the graph in terms of node degree");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> og_file(parser, "FILE", "describe node degree in this graph", {'i', "input"});
    args::Flag summarize(parser, "summarize", "summarize the degree with aggregate statistics", {'S', "summarize"});
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

    if (!og_file) {
        std::cerr << "[odgi::degree] error: please specify a target graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    odgi::graph_t graph;
    assert(argc > 0);
    if (!args::get(og_file).empty()) {
        std::string infile = args::get(og_file);
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
            ifstream f(infile.c_str());
            graph.deserialize(f);
            f.close();
        }
    }

    bool in_parallel = false;
    if (args::get(threads)) {
        in_parallel = true;
        omp_set_num_threads(args::get(threads));
    } else {
        omp_set_num_threads(1);
    }

    if (!summarize) {
        std::cout << "#node.id\tnode.degree" << std::endl;
        graph.for_each_handle(
            [&](const handle_t& handle) {
#pragma omp critical (cout)
                std::cout << graph.get_id(handle) << "\t"
                          << graph.get_degree(handle, false) + graph.get_degree(handle, true)
                          << std::endl;
            }, in_parallel);
    } else {
        uint64_t total_edges = 0;
        uint64_t min_degree = std::numeric_limits<uint64_t>::max();
        uint64_t max_degree = std::numeric_limits<uint64_t>::min();
        graph.for_each_handle(
            [&](const handle_t& handle) {
                uint64_t degree = graph.get_degree(handle, false) + graph.get_degree(handle, true);
                total_edges += degree;
                min_degree = std::min(min_degree, degree);
                max_degree = std::max(max_degree, degree);
            });
        std::cout << "#node.count\tedge.count\tavg.degree\tmin.degree\tmax.degree" << std::endl
                  << graph.get_node_count() << "\t"
                  << total_edges / 2 << "\t" // we double-count edges
                  << (double) total_edges / (double)graph.get_node_count() << "\t"
                  << min_degree << "\t"
                  << max_degree
                  << std::endl;
    }

    return 0;
}

static Subcommand odgi_degree("degree", "node degree information",
                              PIPELINE, 3, main_degree);


}
