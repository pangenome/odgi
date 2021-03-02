#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <queue>
#include <weakly_connected_components.hpp>
#include "handlegraph/path_position_handle_graph.hpp"
namespace odgi {

    using namespace odgi::subcommand;


/// We can accumulate a subgraph without accumulating all the edges between its nodes
/// this helper ensures that we get the full set
    void add_connecting_edges_to_subgraph(const graph_t& source, graph_t& subgraph) {
        subgraph.for_each_handle([&](const handle_t& handle) {
            nid_t id = subgraph.get_id(handle);
            handle_t source_handle = source.get_handle(id, subgraph.get_is_reverse(handle));
            source.follow_edges(source_handle, false, [&](const handle_t& next) {
                nid_t next_id = source.get_id(next);
                if (subgraph.has_node(next_id)) {
                    handle_t subgraph_next = subgraph.get_handle(next_id, source.get_is_reverse(next));
                    if (!subgraph.has_edge(handle, subgraph_next)) {
                        subgraph.create_edge(handle, subgraph_next);
                    }
                }
            });
            source.follow_edges(source_handle, true, [&](const handle_t& prev) {
                nid_t prev_id = source.get_id(prev);
                if (subgraph.has_node(prev_id)) {
                    handle_t subgraph_prev = subgraph.get_handle(prev_id, source.get_is_reverse(prev));
                    if (!subgraph.has_edge(subgraph_prev, handle)) {
                        subgraph.create_edge(subgraph_prev, handle);
                    }
                }
            });
        });
    }


    void algorithms_expand_subgraph_by_steps(const graph_t& source, graph_t& subgraph, const uint64_t& steps, bool forward_only) {
        std::vector<handle_t> curr_handles;
        subgraph.for_each_handle([&](const handle_t& h) {
            curr_handles.push_back(h);
        });
        for (uint64_t i = 0; i < steps && !curr_handles.empty(); ++i) {
            std::vector<handle_t> next_handles;
            for (auto& h : curr_handles) {
                handle_t old_h = source.get_handle(subgraph.get_id(h));
                source.follow_edges(old_h, false, [&](const handle_t& c) {
                    handle_t x;
                    if (!subgraph.has_node(source.get_id(c))) {
                        x = subgraph.create_handle(source.get_sequence(source.get_is_reverse(c)?source.flip(c):c), source.get_id(c));
                        next_handles.push_back(x);
                    } else {
                        x = subgraph.get_handle(source.get_id(c));
                    }
                    if (source.get_is_reverse(c)) {
                        x = subgraph.flip(x);
                    }
                    if (!subgraph.has_edge(h, x)) {
                        subgraph.create_edge(h, x);
                    }
                });
                if (!forward_only) {
                    source.follow_edges(old_h, true, [&](const handle_t& c) {
                        handle_t x;
                        if (!subgraph.has_node(source.get_id(c))) {
                            x = subgraph.create_handle(source.get_sequence(source.get_is_reverse(c)?source.flip(c):c), source.get_id(c));
                            next_handles.push_back(x);
                        } else {
                            x = subgraph.get_handle(source.get_id(c));
                        }
                        if (source.get_is_reverse(c)) {
                            x = subgraph.flip(x);
                        }
                        if (!subgraph.has_edge(x, h)) {
                            subgraph.create_edge(x, h);
                        }
                    });
                }
            }
            curr_handles = std::move(next_handles);
        }
        add_connecting_edges_to_subgraph(source, subgraph);
    }

// Create a subpath name
    string Paths_make_subpath_name(const string& path_name, size_t offset, size_t end_offset) {
        string out_name = path_name + "[" + std::to_string(offset);
        if (end_offset > 0) {
            out_name += "-" + std::to_string(end_offset);
        }
        out_name += "]";
        return out_name;
    }
/// add subpaths to the subgraph, providing a concatenation of subpaths that are discontiguous over the subgraph
/// based on their order in the path position index provided by the source graph
/// will clear any path found in both graphs before writing the new steps into it
/// if subpath_naming is true, a suffix will be added to each path in the subgraph denoting its offset
/// in the source graph (unless the subpath was not cut up at all)
    void algorithms_add_subpaths_to_subgraph(const graph_t& source, graph_t& subgraph,
                                  bool subpath_naming) {
        std::unordered_map<std::string, std::map<uint64_t, handle_t> > subpaths;
        subgraph.for_each_handle([&](const handle_t& h) {
            handlegraph::nid_t id = subgraph.get_id(h);
            if (source.has_node(id)) {
                handle_t handle = source.get_handle(id);
                uint64_t pos = 0;
                source.for_each_step_on_handle(handle, [&](const step_handle_t& step) {
                    path_handle_t path = source.get_path_handle_of_step(step);
                    bool is_rev = source.get_is_reverse(source.get_handle_of_step(step));
                    std::string path_name = source.get_path_name(path);
                    subpaths[path_name][pos] = is_rev ? subgraph.flip(h) : h;
                    pos += source.get_length(handle);
                    return true;
                });
            }
        });

        function<path_handle_t(const string&, bool, size_t)> new_subpath =
                [&subgraph](const string& path_name, bool is_circular, size_t subpath_offset) {
                    string subpath_name = Paths_make_subpath_name(path_name, subpath_offset, 0);
                    if (subgraph.has_path(subpath_name)) {
                        subgraph.destroy_path(subgraph.get_path_handle(subpath_name));
                    }
                    return subgraph.create_path_handle(subpath_name, is_circular);
                };

        for (auto& subpath : subpaths) {
            const std::string& path_name = subpath.first;
            path_handle_t source_path_handle = source.get_path_handle(path_name);
            // destroy the path if it exists
            if (subgraph.has_path(path_name)) {
                subgraph.destroy_path(subgraph.get_path_handle(path_name));
            }
            // create a new path.  give it a subpath name if the flag's on and its smaller than original
            path_handle_t path;
            if (!subpath_naming || subpath.second.size() == source.get_step_count(source_path_handle) ||
                subpath.second.empty()) {
                path = subgraph.create_path_handle(path_name, source.get_is_circular(source_path_handle));
            } else {
                path = new_subpath(path_name, source.get_is_circular(source_path_handle), subpath.second.begin()->first);
            }
            for (auto p = subpath.second.begin(); p != subpath.second.end(); ++p) {
                const handle_t& handle = p->second;
                if (p != subpath.second.begin() && subpath_naming) {
                    auto prev = p;
                    --prev;
                    const handle_t& prev_handle = prev->second;
                    // distance from map
                    size_t delta = p->first - prev->first;
                    // what the distance should be if they're contiguous depends on relative orienations
                    size_t cont_delta = subgraph.get_length(prev_handle);
                    if (delta != cont_delta) {
                        // we have a discontinuity!  we'll make a new path can continue from there
                        assert(subgraph.get_step_count(path) > 0);
                        path = new_subpath(path_name, subgraph.get_is_circular(path), p->first);
                    }
                }
                //fill in the path information
                subgraph.append_step(path, handle);
            }
        }
    }



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
        args::Flag debug(parser, "debug", "print information about the components and the progress to stderr",
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

        uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

        std::string output_dir = ".";

        std::vector<ska::flat_hash_set<handlegraph::nid_t>> weak_components = algorithms::weakly_connected_components(&graph);

        for (uint64_t component_index; component_index < weak_components.size(); ++component_index) {
            auto& weak_component = weak_components[component_index];

            graph_t subgraph;

            for (auto node_id : weak_component) {
                subgraph.create_handle(graph.get_sequence(graph.get_handle(node_id)), node_id);
            }

            algorithms_expand_subgraph_by_steps(graph, subgraph, numeric_limits<uint64_t>::max(), false);
            algorithms_add_subpaths_to_subgraph(graph, subgraph, false);

            // Save the component
            string filename = output_dir + "/component" + to_string(component_index) + ".og";

            //todo Now report what paths went into the component in parseable TSV

            ofstream f(filename);
            subgraph.serialize(f);
            f.close();

            std::cerr << component_index << ": " << subgraph.get_node_count() << " - " << subgraph.get_path_count() << std::endl;
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
