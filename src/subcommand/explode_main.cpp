#include "subcommand.hpp"

#include "args.hxx"
#include <queue>
#include <atomic_bitvector.hpp>
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
                "Breaks a graph into connected components storing each component in its own file.");
        args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
        args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
        args::Group explode_opts(parser, "[ Explode Options ]");
        args::Flag _to_gfa(explode_opts, "to_gfa", "Write each connected component to a file in GFAv1 format.", {'g', "to-gfa"});
        args::ValueFlag<std::string> _prefix(explode_opts, "STRING",
                                             "Write each connected component to a file with the given STRING prefix. "
                                             "The file for the component number `i` will be named `STRING.i.EXTENSION` "
                                             "(default: `component.i.og` or `component.i.gfa`).", {'p', "prefix"});
        args::ValueFlag<uint64_t> _write_biggest_components(explode_opts, "N",
                                                            "Specify the number of the biggest connected components to write, sorted by decreasing size (default: disabled, for writing them all).",
                                                            {'b', "biggest"});
        args::ValueFlag<char> _size_metric(explode_opts, "C",
                                           "Specify how to sort the connected components by size:\np) Path mass (total number of path bases) (default).\nl) Graph length (number of node bases).\nn) Number of nodes.\nP) Longest path.",
                                           {'s', "sorting-criteria"});
        args::Flag _optimize(explode_opts, "optimize", "Compact the node ID space in each connected component.",
                             {'O', "optimize"});
        args::Group threading_opts(parser, "[ Threading ]");
        args::ValueFlag<uint64_t> nthreads(threading_opts, "N",
                                           "Number of threads to use for parallel operations.",
                                           {'t', "threads"});
        args::Group processing_info_opts(parser, "[ Processing Information ]");
        args::Flag _progress(processing_info_opts, "progress", "Print information about the components and the progress to stderr.",
                          {'P', "progress"});
        args::Group program_info_opts(parser, "[ Program Information ]");
        args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi explode.", {'h', "help"});

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

        const bool to_gfa = args::get(_to_gfa);
        const bool optimize = args::get(_optimize);
        const bool progress = args::get(_progress);
        const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

		graph_t graph;
        assert(argc > 0);
        if (!args::get(dg_in_file).empty()) {
            std::string infile = args::get(dg_in_file);
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
                utils::handle_gfa_odgi_input(infile, "explode", progress, num_threads, graph);
            }
        }

        std::string output_dir_plus_prefix = "./";
        if (!args::get(_prefix).empty()) {
            output_dir_plus_prefix += args::get(_prefix);
        } else {
            output_dir_plus_prefix += "component";
        }

        std::vector<ska::flat_hash_set<handlegraph::nid_t>> weak_components =
                algorithms::weakly_connected_components(&graph);


        atomicbitvector::atomic_bv_t ignore_component(weak_components.size());

        if (_write_biggest_components && args::get(_write_biggest_components) > 0) {
            char size_metric = _size_metric ? args::get(_size_metric) : 'p';

            auto get_path_handles = [](const graph_t &graph, const ska::flat_hash_set<handlegraph::nid_t> &node_ids,
                                       set<path_handle_t> &paths) {
                for (auto node_id : node_ids) {
                    handle_t handle = graph.get_handle(node_id);

                    graph.for_each_step_on_handle(handle, [&](const step_handle_t &source_step) {
                        paths.insert(graph.get_path_handle_of_step(source_step));
                    });
                }
            };

            auto get_path_length = [](const graph_t &graph, const path_handle_t &path_handle) {
                uint64_t path_len = 0;
                graph.for_each_step_in_path(path_handle, [&](const step_handle_t &s) {
                    path_len += graph.get_length(graph.get_handle_of_step(s));
                });
                return path_len;
            };

            std::vector<std::pair<uint64_t, uint64_t>> component_and_size;
            component_and_size.resize(weak_components.size());

            // Fill the vector with component sizes
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
            for (uint64_t component_index = 0; component_index < weak_components.size(); ++component_index) {
                ignore_component.set(component_index);

                component_and_size[component_index].first = component_index;

                auto &weak_component = weak_components[component_index];

                uint64_t size = 0;

                switch (size_metric) {
                    case 'l': {
                        // graph length (number of node bases)

                        for (auto node_id : weak_component) {
                            size += graph.get_length(graph.get_handle(node_id));
                        }

                        break;
                    }
                    case 'n': {
                        // number of nodes

                        size = weak_component.size();

                        break;
                    }
                    case 'P': {
                        // longest path",

                        set<path_handle_t> paths;
                        get_path_handles(graph, weak_component, paths);

                        uint64_t current_path_len;
                        for (path_handle_t path_handle : paths) {
                            current_path_len = get_path_length(graph, path_handle);
                            if (current_path_len > size) {
                                size = current_path_len;
                            }
                        }

                        break;
                    }
                    default: {
                        // p: path mass (total number of path bases)

                        set<path_handle_t> paths;
                        get_path_handles(graph, weak_component, paths);

                        for (path_handle_t path_handle : paths) {
                            size += get_path_length(graph, path_handle);
                        }

                        break;
                    }
                }

                component_and_size[component_index].second = size;
            }

            // Sort by component size
            std::sort(component_and_size.begin(), component_and_size.end(), [](auto &a, auto &b) {
                return a.second > b.second;
            });

            // Not ignore the first `write_biggest_components` components
            uint64_t write_biggest_components = args::get(_write_biggest_components);
            for (auto &c_and_s : component_and_size) {
                ignore_component.reset(c_and_s.first);

                if (--write_biggest_components == 0) {
                    break;
                }
            }

//            for(auto& c : component_and_size) {
//                std::cerr << c.first << " (" << ignore_component.test(c.first) << ") - " << c.second << std::endl;
//            }
        }

        std::unique_ptr<algorithms::progress_meter::ProgressMeter> component_progress;
        if (progress) {
            component_progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                    weak_components.size(), "[odgi::explode] exploding component(s)");

            std::cerr << "[odgi::explode] detected " << weak_components.size() << " connected component(s)"
                      << std::endl;

            if (_write_biggest_components && args::get(_write_biggest_components) > 0) {
                uint64_t write_biggest_components = args::get(_write_biggest_components);

                std::cerr << "[odgi::explode] explode the "
                          << (write_biggest_components <= weak_components.size() ? write_biggest_components
                                                                                 : weak_components.size())
                          << " biggest connected component(s)" << std::endl;
            }
        }

        for (uint64_t component_index = 0; component_index < weak_components.size(); ++component_index) {
            if (!ignore_component.test(component_index)) {
                auto &weak_component = weak_components[component_index];

                graph_t subgraph;

                for (auto node_id : weak_component) {
                    subgraph.create_handle(graph.get_sequence(graph.get_handle(node_id)), node_id);
                }

                weak_component.clear();

                algorithms::add_connecting_edges_to_subgraph(graph, subgraph);
                algorithms::add_full_paths_to_component(graph, subgraph, num_threads);

                if (optimize) {
                    subgraph.optimize();
                }

                const string filename = output_dir_plus_prefix + "." + to_string(component_index) + (to_gfa ? ".gfa" : ".og");

                // Save the component
                ofstream f(filename);
                if (to_gfa){
                    subgraph.to_gfa(f, false);
                }else {
                    subgraph.serialize(f);
                }
                f.close();

                /*if (debug) {
#pragma omp critical (cout)
                        std::cerr << "Written component num. " << component_index
                                  << " - num. of nodes " << subgraph.get_node_count()
                                  << " - num. of paths: " << subgraph.get_path_count()
                                  << std::endl;
                }*/
            }

            if (progress) {
                component_progress->increment(1);
            }
        }

        if (progress) {
            component_progress->finish();
        }

        return 0;
    }

    static Subcommand odgi_explode("explode", "Breaks a graph into connected components storing each component in its own file.",
                                   PIPELINE, 3, main_explode);

}
