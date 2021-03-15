#include "subcommand.hpp"

#include "args.hxx"
#include <queue>

#include "src/algorithms/subgraph/extract.hpp"

namespace odgi {

    using namespace odgi::subcommand;

    int main_explode(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi explode";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser(
                "breaks a graph into connected components in their own files");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
        args::ValueFlag<std::string> _prefix(parser, "STRING",
                                             "write each connected component in a file with the given prefix. "
                                             "The file for the component `i` will be named `STRING.i.og` "
                                             "(default: `component`)\"", {'p', "prefix"});
        args::Flag _optimize(parser, "optimize", "compact the node ID space in each connected component",
                             {'O', "optimize"});
        args::ValueFlag<uint64_t> nthreads(parser, "N", "number of threads to use (to write the components in parallel)", {'t', "threads"});
        args::Flag _debug(parser, "progress", "print information about the components and the progress to stderr",
                          {'P', "progress"});

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
                    << "[odgi::explode] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
                    << std::endl;
            return 1;
        }

        graph_t graph;
        assert(argc > 0);
        std::string infile = args::get(dg_in_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
                ifstream f(infile.c_str());
                graph.deserialize(f);
                f.close();
            }
        }

        bool debug = args::get(_debug);
        bool optimize = args::get(_optimize);

        uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

        std::string output_dir_plus_prefix = "./";
        if (!args::get(_prefix).empty()) {
            output_dir_plus_prefix += args::get(_prefix);
        } else {
            output_dir_plus_prefix += "component";
        }

        std::vector<ska::flat_hash_set<handlegraph::nid_t>> weak_components =
                algorithms::weakly_connected_components(&graph);

        std::unique_ptr<algorithms::progress_meter::ProgressMeter> component_progress;
        if (debug) {
            component_progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                    weak_components.size(), "[odgi::explode] exploding components");

            std::cerr << "[odgi::explode] detected " << weak_components.size() << " connected components" << std::endl;
        }

        //std::mutex debug_mutex;

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
        for (uint64_t component_index = 0; component_index < weak_components.size(); ++component_index) {
            auto &weak_component = weak_components[component_index];

            graph_t subgraph;

            for (auto node_id : weak_component) {
                subgraph.create_handle(graph.get_sequence(graph.get_handle(node_id)), node_id);
            }

            algorithms::add_connecting_edges_to_subgraph(graph, subgraph);
            algorithms::add_full_paths_to_component(graph, subgraph);

            if (optimize) {
                subgraph.optimize();
            }

            string filename = output_dir_plus_prefix + "." + to_string(component_index) + ".og";

            // Save the component
            ofstream f(filename);
            subgraph.serialize(f);
            f.close();

            /*if (debug) {
                {
                    std::lock_guard<std::mutex> guard(debug_mutex);

                    std::cerr << "Written component num. " << component_index
                              << " - num. of nodes " << subgraph.get_node_count()
                              << " - num. of paths: " << subgraph.get_path_count()
                              << std::endl;
                }
            }*/

            if (debug) {
                component_progress->increment(1);
            }
        }

        if (debug) {
            component_progress->finish();
        }

        return 0;
    }

    static Subcommand odgi_explode("explode", "breaks a graph into connected components",
                                   PIPELINE, 3, main_explode);


}
