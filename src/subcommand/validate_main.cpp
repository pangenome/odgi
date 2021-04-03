#include "subcommand.hpp"
#include "odgi.hpp"
#include "position.hpp"
#include "args.hxx"
#include "split.hpp"
#include "algorithms/bfs.hpp"
#include <omp.h>

namespace odgi {

    using namespace odgi::subcommand;

    int main_validate(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi validate";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser(
                "validate the graph (currently, check if the paths are consistent with the graph topology)");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> og_file(parser, "FILE", "validate this graph", {'i', "input"});

        args::ValueFlag<uint64_t> nthreads(parser, "N", "number of threads to use", {'t', "threads"});

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

        if (!og_file) {
            std::cerr << "[odgi::validate] error: please specify a graph to validate via -i=[FILE], --idx=[FILE]."
                      << std::endl;
            return 1;
        }

        odgi::graph_t graph;
        assert(argc > 0);
        std::string infile = args::get(og_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
                ifstream f(infile.c_str());
                graph.deserialize(f);
                f.close();
            }
        }

        uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;
        omp_set_num_threads(num_threads);

        bool valid_graph = true;

        std::vector<path_handle_t> paths;
        paths.reserve(graph.get_path_count());
        graph.for_each_path_handle([&](const path_handle_t &path) {
            paths.push_back(path);
        });

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
        for (auto path : paths) {
            graph.for_each_step_in_path(path, [&](const step_handle_t &step) {
                if (graph.has_next_step(step)) {
                    step_handle_t next_step = graph.get_next_step(step);
                    handle_t h = graph.get_handle_of_step(step);
                    handle_t next_h = graph.get_handle_of_step(next_step);

                    if (!graph.has_edge(h, next_h)) {
#pragma omp critical (cout)
                        std::cerr << "[odgi::validate] error: the path " << graph.get_path_name(path) << " does not "
                                  << "respect the graph topology: the link "
                                  << graph.get_id(h) << (graph.get_is_reverse(h) ? "-" : "+")
                                  << ","
                                  << graph.get_id(next_h) << (graph.get_is_reverse(next_h) ? "-" : "+")
                                  << " is missing." << std::endl;

                        valid_graph = false;
                    }
                }
            });
        }

        return (valid_graph ? 0 : 1);
    }

    static Subcommand odgi_validate("validate",
                                    "validate the graph (currently, check if the paths are consistent with the graph topology)",
                                    PIPELINE, 3, main_validate);

}
