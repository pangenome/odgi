#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/cover.hpp"

namespace odgi {

    using namespace odgi::subcommand;

    int main_cover(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi cover";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser("find a path cover of the graph");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
        args::ValueFlag<std::string> dg_out_file(parser, "FILE",
                                                 "store the graph with the generated paths in this file", {'o', "out"});
        args::ValueFlag<uint64_t> num_paths_per_component(parser, "N", "number of paths to generate per component",
                                                          {'n', "num-paths-per-component"});
        args::ValueFlag<uint64_t> node_window_size(parser, "N",
                                                   "size of the node window to check each time a new path is extended (it has to be greater than or equal to 2)",
                                                   {'k', "node-window-size"});
        args::ValueFlag<uint64_t> min_node_coverage(parser, "N",
                                                    "minimum node coverage to reach (it has to be greater than 0)",
                                                    {'c', "min-node-coverage"});
        args::Flag debug(parser, "debug", "print information about the components and the progress to stdout",
                         {'d', "debug"});

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

        if (!dg_in_file) {
            std::cerr
                    << "[odgi cover] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
                    << std::endl;
            return 1;
        }

        if (!dg_out_file) {
            std::cerr
                    << "[odgi cover] error: please specify an output file to where to store the graph via -o=[FILE], --out=[FILE]."
                    << std::endl;
            return 1;
        }

        uint64_t _num_paths_per_component = args::get(num_paths_per_component);
        uint64_t _min_node_coverage = args::get(min_node_coverage);
        if (_num_paths_per_component && _min_node_coverage) {
            std::cerr
                    << "[odgi cover] error: please specify -n/--num-paths-per-component or -c/--min-node-coverage, not both."
                    << std::endl;
            return 1;
        } else if (!_num_paths_per_component && !_min_node_coverage) {
            _num_paths_per_component = algorithms::PATH_COVER_DEFAULT_N;
        }

        uint64_t _node_window_size = args::get(node_window_size) ? args::get(node_window_size)
                                                                 : algorithms::PATH_COVER_DEFAULT_K;
        if (_node_window_size < 2) {
            std::cerr
                    << "[odgi cover] error: please specify a node window size greater than or equal to 2 via -k=[N], --node-window-size=[N]."
                    << std::endl;
            return 1;
        }

        graph_t graph;
        assert(argc > 0);
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

        uint64_t max_number_of_paths_generable = graph.get_node_count();
        if (_min_node_coverage) {
            std::cerr << "There will be generated paths until the minimum node coverage is " << _min_node_coverage
                      << ", or until the maximum number of allowed generated paths is reached ("
                      << max_number_of_paths_generable << ")." << std::endl;
        } else {
            std::cerr << "There will be generated " << _num_paths_per_component << " paths per component." << std::endl;
        }

        algorithms::path_cover(graph, _num_paths_per_component, _node_window_size, _min_node_coverage,
                               max_number_of_paths_generable, args::get(debug));

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

    static Subcommand odgi_cover("cover", "find a path cover of the graph",
                                 PIPELINE, 3, main_cover);


}
