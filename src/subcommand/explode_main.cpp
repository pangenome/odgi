#include "subcommand.hpp"

#include "args.hxx"
#include <queue>

#include "algorithms/explode.hpp"

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
                "breaks a graph into connected components in their own files in the given directory");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
        args::ValueFlag<std::string> dg_out_file(parser, "FILE",
                                                 "store the graph with the generated paths in this file", {'o', "out"});
        args::ValueFlag<uint64_t> nthreads(parser, "N", "number of threads to use", {'t', "threads"});
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

        if (!dg_in_file) {
            std::cerr
                    << "[odgi::explode] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
                    << std::endl;
            return 1;
        }

        if (!dg_out_file) {
            std::cerr
                    << "[odgi::explode] error: please specify an output file to where to store the graph via -o=[FILE], --out=[FILE]."
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

        bool debug = args::get(_debug);

        uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

        std::string output_dir = ".";

        std::vector<ska::flat_hash_set<handlegraph::nid_t>> weak_components = algorithms::weakly_connected_components(
                &graph);

        for (uint64_t component_index; component_index < weak_components.size(); ++component_index) {
            auto &weak_component = weak_components[component_index];

            graph_t subgraph;

            for (auto node_id : weak_component) {
                subgraph.create_handle(graph.get_sequence(graph.get_handle(node_id)), node_id);
            }

            algorithms::expand_subgraph_by_steps(graph, subgraph, numeric_limits<uint64_t>::max(), false);
            algorithms::add_subpaths_to_subgraph(graph, subgraph, false);

            // Save the component
            string filename = output_dir + "/component" + to_string(component_index) + ".og";

            //todo Now report what paths went into the component in parseable TSV

            ofstream f(filename);
            subgraph.serialize(f);
            f.close();

            if (debug) {
                std::cerr << "Written component num. " << component_index
                          << " - num. of nodes " << subgraph.get_node_count()
                          << " - num. of paths: " << subgraph.get_path_count()
                          << std::endl;
            }

        }






        /*
        // Count through the components we build
        size_t component_index = 0;

        // Track all the nodes we've already assigned to subgraphs
        unordered_set<handle_t> node_assigned;

        graph.for_each_handle([&](const handle_t &start) {
            if (!node_assigned.count(start)) {
                std::cerr << "component_index " << component_index << std::endl;

                // It's a new connected component!
                graph_t component;

                // We want to track the path names in each component
                set<path_handle_t> paths;

                deque<handle_t> queue{start};

                // Mark this connected node as used in a component.
                node_assigned.insert(start);

                while (!queue.empty()) {
                    handle_t handle = queue.front();
                    queue.pop_front();

                    // Copy node over
                    handle_t new_handle = component.create_handle(graph.get_sequence(handle), graph.get_id(handle));

                    //std::cerr << "A " << graph.get_id(handle) << std::endl;
                    //std::cerr << "B " << component.get_id(new_handle) << std::endl;

                    // Copy over its edges and queue the next handles
                    graph.follow_edges(handle, false, [&](const handle_t &next) {
                        //std::cerr << "next " << component.get_id(next) << std::endl;

                        if (component.has_node(graph.get_id(next))) {
                            component.create_edge(new_handle,
                                                  component.get_handle(graph.get_id(next), graph.get_is_reverse(next)));
                        }
                        if (!node_assigned.count(next)) {
                            queue.push_back(next);
                            node_assigned.insert(next);
                        }
                    });
                    graph.follow_edges(handle, true, [&](const handle_t &prev) {
                        //std::cerr << "prev " << component.get_id(prev) << std::endl;

                        if (component.has_node(graph.get_id(prev))) {
                            component.create_edge(component.get_handle(graph.get_id(prev),
                                                                       graph.get_is_reverse(prev)), new_handle);
                        }
                        if (!node_assigned.count(prev)) {
                            queue.push_back(prev);
                            node_assigned.insert(prev);
                        }
                    });

                    // Record paths
                    graph.for_each_step_on_handle(handle, [&](const step_handle_t &step) {
                        paths.insert(graph.get_path_handle_of_step(step));
                    });
                }

                //std::cerr << "C" << std::endl;

                // Copy the paths over
                for (path_handle_t path_handle : paths) {
                    path_handle_t new_path_handle = component.create_path_handle(graph.get_path_name(path_handle),
                                                                                 graph.get_is_circular(path_handle));
                    for (handle_t handle : graph.scan_path(path_handle)) {
                        component.append_step(new_path_handle, component.get_handle(graph.get_id(handle),
                                                                                    graph.get_is_reverse(handle)));
                    }
                }

                // Save the component
                string filename = output_dir + "/component" + to_string(component_index) + ".og";

                // Now report what paths went into the component in parseable TSV
                cout << filename;
                for (auto &path_handle : paths) {
                    cout << "\t" << graph.get_path_name(path_handle);
                }
                cout << endl;

                ofstream f(filename);
                component.serialize(f);
                f.close();

                component_index++;
            }
        });
*/

//        std::string outfile = args::get(dg_out_file);
//        if (outfile.size()) {
//            if (outfile == "-") {
//                graph.serialize(std::cout);
//            } else {
//                ofstream f(outfile.c_str());
//                graph.serialize(f);
//                f.close();
//            }
//        }
        return 0;
    }

    static Subcommand odgi_explode("explode", "breaks a graph into connected components",
                                   PIPELINE, 3, main_explode);


}
