#include "subcommand.hpp"

#include "odgi.hpp"
#include "handlegraph/path_position_handle_graph.hpp"
#include "progress.hpp"

#include "args.hxx"
#include <queue>

#include "utils.hpp"

namespace odgi {

    using namespace odgi::subcommand;

    int main_squeeze(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi squeeze";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser(
                "Squeezes multiple graphs in ODGI format into the same file in ODGI format.");
        args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
        args::ValueFlag<std::string> _input_graphs(mandatory_opts, "FILE",
                                                   "Input file containing the list of graphs to squeeze into the same\n"
                                                   "  file. The file must contain one graph per line. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!",
                                                   {'f', "input-graphs"});
        args::ValueFlag<std::string> dg_out_file(mandatory_opts, "FILE", "Store all the input graphs in this file. The file name usually ends with *.og*.",
                                                 {'o', "out"});
        args::Group squeeze_opt(parser, "[ Squeeze Options ]");
        args::ValueFlag<char> _add_suffix(squeeze_opt, "STRING",
                                          "Add the separator and the input file rank as suffix to the path names\n"
                                          "  (to avoid path name collisions).",
                                          {'s', "rank-suffix"});
        args::Flag _optimize(parser, "optimize", "Compact the node ID space for each connected component before squeezing.",
                             {'O', "optimize"});
        args::Group threading_opts(parser, "[ Threading ]");
        args::ValueFlag<uint64_t> nthreads(threading_opts, "N", "Number of threads to use for parallel operations.",
                                           {'t', "threads"});
        args::Group processing_info_opts(parser, "[ Processing Information ]");
        args::Flag progress(parser, "progress", "Print information about the progress to stderr.",
                          {'P', "progress"});
        args::Group program_info_opts(parser, "[ Program Information ]");
        args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi squeeze.", {'h', "help"});

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
                    << "[odgi::squeeze] error: please specify an input file from where to take the path of the graphs to squeeze via -i=[FILE], --idx=[FILE]."
                    << std::endl;
            return 1;
        }

        if (!dg_out_file) {
            std::cerr
                    << "[odgi::squeeze] error: please specify an output file to where to store the graphs via -o=[FILE], --out=[FILE]."
                    << std::endl;
            return 1;
        }

        uint64_t num_input_graphs = 0;
        {

            std::ifstream file_input_graphs(input_graphs);
            std::string line;
            while (std::getline(file_input_graphs, line)) {
                if (!line.empty()) {
                    ++num_input_graphs;
                }
            }
            file_input_graphs.close();
        }

        if (num_input_graphs == 0) {
            std::cerr
                    << "[odgi::squeeze] error: the input file contains no input graph paths."
                    << std::endl;
            return 1;
        }

        bool debug = args::get(progress);
        bool optimize = args::get(_optimize);

        uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;
        omp_set_num_threads(num_threads);

        std::unique_ptr<algorithms::progress_meter::ProgressMeter> component_progress;

        char separator;
        if (_add_suffix) {
            separator = args::get(_add_suffix);
        }

        std::unique_ptr<algorithms::progress_meter::ProgressMeter> squeeze_progress;
        if (debug) {
            squeeze_progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                    num_input_graphs, "[odgi::squeeze] squeezing input graphs");

            std::cerr << "[odgi::squeeze] detected " << num_input_graphs << " input graphs" << std::endl;
        }


        uint64_t shift_id = 0;
        graph_t squeezed_graph;

        uint64_t input_graph_rank = 0;
        std::ifstream file_input_graphs(input_graphs);
        std::string line;
        while (std::getline(file_input_graphs, line)) {
            if (!line.empty()) {
                graph_t graph;

				utils::handle_gfa_odgi_input(line, "squeeze", args::get(progress), num_threads, graph);

                if (optimize) {
                    graph.optimize();
                }

                uint64_t max_id = 0;
                uint64_t new_node_id;

                graph.for_each_handle([&](const handle_t &h) {
                    new_node_id = graph.get_id(h) + shift_id;

                    squeezed_graph.create_handle(graph.get_sequence(h), new_node_id);

                    if (new_node_id > max_id) {
                        max_id = new_node_id;
                    }
                });

                // add contacts for the edges
                graph.for_each_handle([&](const handle_t &h) {
                    handle_t new_handle_h = squeezed_graph.get_handle(graph.get_id(h) + shift_id);
                    const bool h_is_rev = graph.get_is_reverse(h);

                    graph.follow_edges(h, false, [&](const handle_t &o) {
                        handle_t new_handle_o = squeezed_graph.get_handle(graph.get_id(o) + shift_id);
                        const bool o_is_rev = graph.get_is_reverse(o);

                        squeezed_graph.create_edge(
                            h_is_rev ? graph.flip(new_handle_h) : new_handle_h,
                            o_is_rev ? graph.flip(new_handle_o) : new_handle_o
                        );
                    });
                    graph.follow_edges(h, true, [&](const handle_t &o) {
                        handle_t new_handle_o = squeezed_graph.get_handle(graph.get_id(o) + shift_id);
                        const bool o_is_rev = graph.get_is_reverse(o);

                        squeezed_graph.create_edge(
                            o_is_rev ? graph.flip(new_handle_o) : new_handle_o,
                            h_is_rev ? graph.flip(new_handle_h) : new_handle_h
                        );
                    });
                });

                // Copy the paths
                std::vector<std::pair<path_handle_t, path_handle_t>> old_and_new_paths;
                old_and_new_paths.reserve(graph.get_path_count());

                graph.for_each_path_handle([&](const path_handle_t old_path_handle) {
                    std::string new_path_name = graph.get_path_name(old_path_handle);
                    if (_add_suffix) {
                        new_path_name += separator + std::to_string(input_graph_rank);
                    }

                    path_handle_t new_path_handle = squeezed_graph.create_path_handle(
                            new_path_name, graph.get_is_circular(old_path_handle));

                    old_and_new_paths.push_back({old_path_handle, new_path_handle});
                });

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
                for (auto old_new_path : old_and_new_paths) {
                    graph.for_each_step_in_path(old_new_path.first, [&](const step_handle_t &step) {
                        handle_t old_handle = graph.get_handle_of_step(step);
                        handle_t new_handle = squeezed_graph.get_handle(
                                graph.get_id(old_handle) + shift_id,
                                graph.get_is_reverse(old_handle));

                        squeezed_graph.append_step(old_new_path.second, new_handle);
                    });
                }

                shift_id = max_id;
                ++input_graph_rank;

                if (debug) {
                    squeeze_progress->increment(1);
                }
            }
        }
        file_input_graphs.close();

        if (debug) {
            squeeze_progress->finish();
        }

        std::string outfile = args::get(dg_out_file);
        if (!outfile.empty()) {
            if (outfile == "-") {
                squeezed_graph.serialize(std::cout);
            } else {
                ofstream f(outfile.c_str());
                squeezed_graph.serialize(f);
                f.close();
            }
        }

        return 0;
    }

    static Subcommand odgi_squeeze("squeeze", "Squeezes multiple graphs in ODGI format into the same file in ODGI format.",
                                   PIPELINE, 3, main_squeeze);


}
