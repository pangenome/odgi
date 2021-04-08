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
                "breaks a graph into connected components in their own files");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
        args::ValueFlag<std::string> _prefix(parser, "STRING",
                                             "write each connected component in a file with the given prefix. "
                                             "The file for the component `i` will be named `STRING.i.og` "
                                             "(default: `component`)", {'p', "prefix"});

        args::ValueFlag<uint64_t> _write_biggest_components(parser, "N",
                                                            "specify the number of the biggest connected components to write, sorted by decreasing size (default disabled, for writing them all) ",
                                                            {'b', "biggest"});
        args::ValueFlag<char> _size_metric(parser, "C",
                                           "specify how to sort the connected components by size:\np) path mass (total number of path bases) (default)\nl) graph length (number of node bases)\nn) number of nodes\nP)longest path",
                                           {'s', "sorting-criteria"});

        args::Flag _optimize(parser, "optimize", "compact the node ID space in each connected component",
                             {'O', "optimize"});
        args::ValueFlag<uint64_t> nthreads(parser, "N",
                                           "number of threads to use (to write the components in parallel)",
                                           {'t', "threads"});
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
        if (!args::get(dg_in_file).empty()) {
            std::string infile = args::get(dg_in_file);
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
                ifstream f(infile.c_str());
                graph.deserialize(f);
                f.close();
            }
        }

        const bool debug = args::get(_debug);
        const bool optimize = args::get(_optimize);

        const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

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
                component_and_size[component_index].second = 0;

                auto &weak_component = weak_components[component_index];

                switch (size_metric) {
                    case 'l': {
                        for (auto node_id : weak_component) {
                            component_and_size[component_index].second += graph.get_length(graph.get_handle(node_id));
                        }
                        break;
                    }
                    case 'n': {
                        component_and_size[component_index].second = weak_component.size();
                        break;
                    }
                    case 'P': {
                        set<path_handle_t> paths;
                        for (auto node_id : weak_component) {
                            handle_t handle = graph.get_handle(node_id);

                            graph.for_each_step_on_handle(handle, [&](const step_handle_t &source_step) {
                                paths.insert(graph.get_path_handle_of_step(source_step));
                            });
                        }

                        uint64_t max_path_len = 0, current_path_len;
                        for (path_handle_t path_handle : paths) {
                            current_path_len = get_path_length(graph, path_handle);
                            if (current_path_len > max_path_len) {
                                max_path_len = current_path_len;
                            }
                        }
                        component_and_size[component_index].second = max_path_len;

                        break;
                    }
                    default: {
                        // p path mass (total number of path bases) (default)
                        //ToDo
                        break;
                    }
                }
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
//            exit(1);
        }

        std::unique_ptr<algorithms::progress_meter::ProgressMeter> component_progress;
        if (debug) {
            component_progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                    weak_components.size(), "[odgi::explode] exploding component(s)");

            std::cerr << "[odgi::explode] detected " << weak_components.size() << " connected component(s)"
                      << std::endl;

            if (_write_biggest_components && args::get(_write_biggest_components) > 0) {
                uint64_t write_biggest_components = args::get(_write_biggest_components);

                std::cerr << "[odgi::explode] explode the biggest "
                          << (write_biggest_components <= weak_components.size() ? write_biggest_components
                                                                                 : weak_components.size())
                          << " connected component(s)" << std::endl;
            }
        }

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
        for (uint64_t component_index = 0; component_index < weak_components.size(); ++component_index) {
            if (!ignore_component.test(component_index)) {
                auto &weak_component = weak_components[component_index];

                graph_t subgraph;

                for (auto node_id : weak_component) {
                    subgraph.create_handle(graph.get_sequence(graph.get_handle(node_id)), node_id);
                }

                weak_component.clear();

                algorithms::add_connecting_edges_to_subgraph(graph, subgraph);
                algorithms::add_full_paths_to_component(graph, subgraph);

                if (optimize) {
                    subgraph.optimize();
                }

                const string filename = output_dir_plus_prefix + "." + to_string(component_index) + ".og";

                // Save the component
                ofstream f(filename);
                subgraph.serialize(f);
                f.close();

                /*if (debug) {
#pragma omp critical (cout)
                        std::cerr << "Written component num. " << component_index
                                  << " - num. of nodes " << subgraph.get_node_count()
                                  << " - num. of paths: " << subgraph.get_path_count()
                                  << std::endl;
                }*/
            }

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
