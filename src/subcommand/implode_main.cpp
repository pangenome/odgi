#include "subcommand.hpp"

#include "args.hxx"
#include <queue>

#include "algorithms/explode.hpp"

namespace odgi {

    using namespace odgi::subcommand;

    int main_implode(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi implode";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser(
                "merges multiple graphs into the same file");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> _input_graphs(parser, "FILE",
                                                   "list of graphs to implode into the same file; the file must contain one path per line.",
                                                   {'i', "graphs-to-implode"});
        args::ValueFlag<char> _add_suffix(parser, "C",
                                          "add the separator and the input file rank as prefix to the path names (to avoid path name collisions)",
                                          {'s', "rank-suffix"});

        args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store all the input graphs in this file",
                                                 {'o', "out"});

        args::Flag _optimize(parser, "optimize", "compact the node ID space in each connected component",
                             {'O', "optimize"});
        args::ValueFlag<uint64_t> nthreads(parser, "N",
                                           "number of threads to use (to write the components in parallel)",
                                           {'t', "threads"});
        args::Flag _debug(parser, "debug", "print information about the components and the progress to stderr",
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

        std::string input_graphs = args::get(_input_graphs);
        if (input_graphs.empty()) {
            std::cerr
                    << "[odgi::implode] error: please specify an input file from where to take the path of the graphs to implode via -i=[FILE], --idx=[FILE]."
                    << std::endl;
            return 1;
        }

        if (!dg_out_file) {
            std::cerr
                    << "[odgi::implode] error: please specify an output file to where to store the graphs via -o=[FILE], --out=[FILE]."
                    << std::endl;
            return 1;
        }

        bool debug = args::get(_debug);
        bool optimize = args::get(_optimize);

        uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

        char separator;
        if (_add_suffix) {
            separator = args::get(_add_suffix);
        }

        uint64_t shift_id = 0;
        graph_t imploded_graph;

        uint64_t input_graph_rank = 0;
        std::ifstream file_input_graphs(input_graphs);
        std::string line;
        while (std::getline(file_input_graphs, line)) {
            if (!line.empty()) {
                graph_t graph;
                assert(argc > 0);

                ifstream f(line.c_str());
                graph.deserialize(f);
                f.close();

                if (optimize) {
                    graph.optimize();
                }

                uint64_t max_id = 0;
                graph.for_each_handle([&](const handle_t &h) {
                    uint64_t node_id = graph.get_id(h);
                    uint64_t new_node_id = node_id + shift_id;
                    max_id = std::max(max_id, new_node_id);

                    imploded_graph.create_handle(graph.get_sequence(graph.get_handle(node_id)), new_node_id);
                });

                handle_t new_handle_h, new_handle_o;

                // add contacts for the edges
                graph.for_each_handle([&](const handle_t &h) {
                    new_handle_h = imploded_graph.get_handle(graph.get_id(h) + shift_id);

                    graph.follow_edges(h, false, [&](const handle_t &o) {
                        new_handle_o = imploded_graph.get_handle(graph.get_id(o) + shift_id);
                        imploded_graph.create_edge(new_handle_h, new_handle_o);
                    });
                    graph.follow_edges(h, true, [&](const handle_t &o) {
                        new_handle_o = imploded_graph.get_handle(graph.get_id(o) + shift_id);
                        imploded_graph.create_edge(new_handle_o, new_handle_h);
                    });
                });

                // Copy the paths over
                graph.for_each_path_handle([&](const path_handle_t &p) {
                    std::string new_path_name = graph.get_path_name(p);
                    if (_add_suffix) {
                        new_path_name += separator + std::to_string(input_graph_rank);
                    }

                    path_handle_t new_path_handle = imploded_graph.create_path_handle(new_path_name,
                                                                                      graph.get_is_circular(p));

                    for (handle_t o : graph.scan_path(p)) {
                        new_handle_o = imploded_graph.get_handle(graph.get_id(o) + shift_id);
                        imploded_graph.append_step(new_path_handle, imploded_graph.get_handle(
                                graph.get_id(new_handle_o),
                                imploded_graph.get_is_reverse(new_handle_o))
                        );
                    }
                });


                shift_id = max_id;
                ++input_graph_rank;
            }
        }
        file_input_graphs.close();

        std::string outfile = args::get(dg_out_file);
        if (!outfile.empty()) {
            if (outfile == "-") {
                imploded_graph.serialize(std::cout);
            } else {
                ofstream f(outfile.c_str());
                imploded_graph.serialize(f);
                f.close();
            }
        }

        return 0;
    }

    static Subcommand odgi_explode("implode", "merges multiple graphs into the same file",
                                   PIPELINE, 3, main_implode);


}
