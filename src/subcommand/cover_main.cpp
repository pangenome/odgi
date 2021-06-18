#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/cover.hpp"
#include "utils.hpp"

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

        args::ArgumentParser parser("Cover the graph with paths.");
        args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
        args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
        args::ValueFlag<std::string> dg_out_file(mandatory_opts, "FILE","Write the succinct variation graph with the generated paths in ODGI format to *FILE*. A file ending with *.og* is recommended.", {'o', "out"});
        args::Group cover_opts(parser, "[ Cover Options ]");
        args::ValueFlag<double> hogwild_depth(cover_opts, "DEPTH", "Randomly cover the graph until it reaches the given average DEPTH. Specifying this option overwrites all other cover options except -I, --ignore-paths!",{'H', "hogwild-depth"});
        args::ValueFlag<uint64_t> num_paths_per_component(cover_opts, "N", "Number of paths to generate per component.",{'n', "num-paths-per-component"});
        args::ValueFlag<uint64_t> node_window_size(cover_opts, "N","Size of the node window to check each time a new path is extended (it has to be greater than or equal to 2).",{'k', "node-window-size"});
        args::ValueFlag<uint64_t> min_node_depth(cover_opts, "N","Minimum node depth to reach (it has to be greater than 0). There will be generated paths until the specified minimum node coverage is reached, or until the maximum number of allowed generated paths is reached (number of nodes in the input ODGI graph).",{'c', "min-node-depth"});
        args::Flag ignore_paths(cover_opts, "ignore-paths", "Ignore the paths already embedded in the graph during the node depth initialization.",{'I', "ignore-paths"});
        args::ValueFlag<std::string> write_node_depth(cover_opts, "FILE","Write the node depth at the end of the paths generation to FILE. The file will contain tab-separated values (header included), and have 3 columns: component_id, node_ide, coverage.",{'w', "write-node-depth"});
        args::Group threading_opts(parser, "[ Threading ]");
        args::ValueFlag<uint64_t> nthreads(threading_opts, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
        args::Group processing_info_opts(parser, "[ Processing Information ]");
        args::Flag debug(processing_info_opts, "progress", "Print information about the components and the progress to stderr",{'P', "progress"});
        args::Group program_info_opts(parser, "[ Program Information ]");
        args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi cover.", {'h', "help"});

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
                    << "[odgi::cover] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
                    << std::endl;
            return 1;
        }

        if (!dg_out_file) {
            std::cerr
                    << "[odgi::cover] error: please specify an output file to where to store the graph via -o=[FILE], --out=[FILE]."
                    << std::endl;
            return 1;
        }

        uint64_t _num_paths_per_component = args::get(num_paths_per_component);
        const uint64_t _min_node_depth = args::get(min_node_depth);
        if (_num_paths_per_component && _min_node_depth) {
            std::cerr
                    << "[odgi::cover] error: please specify -n/--num-paths-per-component or -c/--min-node-depth, not both."
                    << std::endl;
            return 1;
        } else if (!_num_paths_per_component && !_min_node_depth) {
            _num_paths_per_component = algorithms::PATH_COVER_DEFAULT_N;
        }

        const uint64_t _node_window_size = args::get(node_window_size) ? args::get(node_window_size)
                                                                 : algorithms::PATH_COVER_DEFAULT_K;
        if (_node_window_size < 2) {
            std::cerr
                    << "[odgi::cover] error: please specify a node window size greater than or equal to 2 via -k=[N], --node-window-size=[N]."
                    << std::endl;
            return 1;
        }

		const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

		graph_t graph;
        assert(argc > 0);
        if (!args::get(dg_in_file).empty()) {
            std::string infile = args::get(dg_in_file);
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
				utils::handle_gfa_odgi_input(infile, "cover", args::get(debug), num_threads, graph);
            }
        }

        if (hogwild_depth) {
            algorithms::hogwild_path_cover(graph, args::get(hogwild_depth), num_threads, args::get(ignore_paths), args::get(debug));
        } else {
            uint64_t max_number_of_paths_generable = graph.get_node_count() * 5;
            if (args::get(debug)){
                if (_min_node_depth) {
                    std::cerr << "[odgi::cover] there will be generated paths until the minimum node depth is " << _min_node_depth
                              << ", or until the maximum number of allowed generated paths is reached ("
                              << max_number_of_paths_generable << ")." << std::endl;
                } else {
                    std::cerr << "[odgi::cover] there will be generated " << _num_paths_per_component << " paths per component."
                              << std::endl;
                }
            }

            std::string node_depth;
            algorithms::path_cover(graph, _num_paths_per_component, _node_window_size, _min_node_depth,
                                   max_number_of_paths_generable,
                                   write_node_depth, node_depth,
                                   num_threads, args::get(ignore_paths), args::get(debug));

            if (write_node_depth) {
                std::string covfile = args::get(write_node_depth);

                ofstream f(covfile.c_str());
                f << node_depth;
                f.close();
            }
        }

        if (!args::get(dg_out_file).empty()) {
            std::string outfile = args::get(dg_out_file);
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

    static Subcommand odgi_cover("cover", "Cover the graph with paths.",
                                 PIPELINE, 3, main_cover);


}
