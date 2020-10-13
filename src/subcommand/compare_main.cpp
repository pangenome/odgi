#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "gfa_to_handle.hpp"

namespace odgi {

    using namespace odgi::subcommand;

    int main_compare(int argc, char **argv) {

        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        const std::string prog_name = "odgi compare";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser("find out if an odgi graph equals a GFA graph");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> odgi_in_file(parser, "FILE", "load the odgi graph from this file", {'i', "idx"});
        args::ValueFlag<std::string> gfa_in_file(parser, "FILE", "store the index in this file", {'g', "gfa"});
        args::Flag progress(parser, "progress", "show progress updates", {'p', "progress"});

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
        if (argc == 1) {
            std::cout << parser;
            return 1;
        }

        if (!odgi_in_file) {
            std::cerr
                    << "[odgi compare::main] error: Please specify an odgi input graph from where to load the graph via -i=[FILE], --idx=[FILE]."
                    << std::endl;
            return 1;
        }

        if (!gfa_in_file) {
            std::cerr
                    << "[odgi compare::main] error: Please specify an GFA file from where to load the graph via -g=[FILE], --gfa=[FILE]."
                    << std::endl;
            return 1;
        }

        // read in the odgi graph
        graph_t odgi_graph;
        assert(argc > 0);
        const std::string odgi_in = args::get(odgi_in_file);
        if (odgi_in.size()) {
            if (odgi_in == "-") {
                odgi_graph.deserialize(std::cin);
            } else {
                ifstream f(odgi_in.c_str());
                odgi_graph.deserialize(f);
                f.close();
            }
        }
        // read in the gfa graph
        graph_t gfa_graph;
        std::string gfa_in = args::get(gfa_in_file);
        if (gfa_in.size()) {
            gfa_to_handle(gfa_in, &gfa_graph, args::get(progress));
        }

        bool equal = odgi_graph.compare(gfa_graph);

        if (equal) {
            std::cerr << "[odgi compare::main] SUCCESS: Both graphs are equal!" << std::endl;
        } else {
            std::cerr << "[odgi compare::main] FAILURE: The entered graphs are not equal!" << std::endl;
        }

        return 0;
    }

    static Subcommand odgi_compare("compare", "create a path index for a given graph",
                                     PIPELINE, 3, main_compare);

}
