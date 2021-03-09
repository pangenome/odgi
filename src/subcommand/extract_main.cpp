#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include <regex>
#include "utils.hpp"
#include "src/algorithms/subgraph/region.hpp"
#include "src/algorithms/explode.hpp"

namespace odgi {

    using namespace odgi::subcommand;

    std::pair<step_handle_t, int64_t>
    get_graph_pos(const odgi::graph_t &graph, const path_handle_t path, const uint64_t &pos) {
        auto path_end = graph.path_end(path);
        uint64_t walked = 0;
        for (step_handle_t s = graph.path_begin(path);
             s != path_end; s = graph.get_next_step(s)) {
            handle_t h = graph.get_handle_of_step(s);
            uint64_t node_length = graph.get_length(h);
            if (walked + node_length > pos) {
                return std::make_pair(s, pos - walked);
            }
            walked += node_length;
        }
#pragma omp critical (cout)
        std::cerr << "[odgi::get_graph_pos] warning: position " << graph.get_path_name(path) << ":" << pos
                  << " outside of path" << std::endl;
        return std::make_pair(path_end, -1);
    }

    // Create a subpath name
    string make_subpath_name(const string& path_name, size_t offset, size_t end_offset) {
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
    void algorithms_add_subpaths_to_subgraph(const graph_t& source, graph_t& subgraph, bool subpath_naming) {
        auto get_position_of_step =
                [](const odgi::graph_t& graph, const step_handle_t& step_to_find) {
                    path_handle_t path = graph.get_path_handle_of_step(step_to_find);
                    auto path_end = graph.path_end(path);
                    uint64_t walked = 0;
                    for (step_handle_t s = graph.path_begin(path); s != path_end; s = graph.get_next_step(s)) {
                        if (s == step_to_find) {
                            return walked;
                        }

                        handle_t h = graph.get_handle_of_step(s);
                        uint64_t node_length = graph.get_length(h);
                        walked += node_length;
                    }
#pragma omp critical (cout)
                    std::cerr << "[odgi::get_position_of_step] warning: step " << graph.get_id((graph.get_handle_of_step(step_to_find))) << " in " << graph.get_path_name(path) << " not found" << std::endl;

                    return walked;
                };

        std::unordered_map<std::string, std::map<uint64_t, handle_t> > subpaths;
        subgraph.for_each_handle([&](const handle_t& h) {
            handlegraph::nid_t id = subgraph.get_id(h);
            if (source.has_node(id)) {
                handle_t handle = source.get_handle(id);
                source.for_each_step_on_handle(handle, [&](const step_handle_t& step) {
                    path_handle_t path = source.get_path_handle_of_step(step);
                    std::string path_name = source.get_path_name(path);
                    uint64_t pos = get_position_of_step(source, step);
                    subpaths[path_name][pos] = source.get_is_reverse(source.get_handle_of_step(step)) ? subgraph.flip(h) : h;

                });
                /*//return iteratee(step, get_is_reverse(get_handle_of_step(step)) != get_is_reverse(handle), get_position_of_step(step));
                source.for_each_step_position_on_handle(handle, [&](const step_handle_t& step, const bool& is_rev, const uint64_t& pos) {
                    path_handle_t path = source.get_path_handle_of_step(step);
                    std::string path_name = source.get_path_name(path);
                    subpaths[path_name][pos] = is_rev ? subgraph.flip(h) : h;
                    return true;
                });*/
            }
        });

        function<path_handle_t(const string&, bool, size_t)> new_subpath =
                [&subgraph](const string& path_name, bool is_circular, size_t subpath_offset) {
                    string subpath_name = make_subpath_name(path_name, subpath_offset, 0);
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

    void algorithms_extract_path_range(const graph_t &source, path_handle_t path_handle, int64_t start, int64_t end,
                                       graph_t &subgraph) {

        // ToDo: could be improved getting and extracting simultaneously (not first get_graph_pos and then extract_path_range)?
        std::pair<step_handle_t, uint64_t> start_step_and_start_pos = get_graph_pos(source, path_handle, start);

        step_handle_t start_step = start_step_and_start_pos.first;  //source.get_step_at_position(path_handle, start);
        size_t start_position = start_step_and_start_pos.second;    //source.get_position_of_step(start_step);
        size_t size_needed = end < 0 ? numeric_limits<size_t>::max() : end - start + 1 + start - start_position;
        size_t running_length = 0;

        auto path_end = source.path_end(path_handle);
        for (step_handle_t cur_step = start_step;
             cur_step != path_end && running_length < size_needed; cur_step = source.get_next_step(cur_step)) {

            handle_t cur_handle = source.get_handle_of_step(cur_step);
            subgraph.create_handle(source.get_sequence(cur_handle), source.get_id(cur_handle));
            if (cur_step != start_step) {
                handle_t prev_handle = source.get_handle_of_step(source.get_previous_step(cur_step));
                subgraph.create_edge(
                        subgraph.get_handle(source.get_id(prev_handle), source.get_is_reverse(prev_handle)),
                        subgraph.get_handle(source.get_id(cur_handle), source.get_is_reverse(cur_handle)));
            }
            running_length += source.get_length(cur_handle);

        }
    }

    int main_extract(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi extract";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        //////////////////////////////////////////////////////////////////////////////
        /*
        1) various ways of collecting a subgraph...  either a node list (-l/--node-list), a set of path ranges, a range
         in the sorted graph, or path ranges projected into the sorted graph (what you're doing)
        2) expand this by some context in terms of bp, steps, etc., if desired, or (achieving your goal) get everything
         in the sort order before the beginning and ending node you've collected
        3) collect the path steps on all nodes in the subgraph, then walk them as far as you can go while staying inside
         the subgraph
        4) remove duplicates
        5) rewrite these subpaths as new paths in the output graph, with the "samtools" or "bedtools" range syntax, we
         add a suffix like :n-m where n= the start and m = the end in 0-based half-open format (bedtools numbers)
        */
        //////////////////////////////////////////////////////////////////////////////

        args::ArgumentParser parser("extract parts of the graph as defined by query criteria");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});

        args::ValueFlag<std::string> og_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
        args::ValueFlag<std::string> og_out_file(parser, "FILE", "store the graph self index in this file",
                                                 {'o', "out"});

        args::ValueFlag<std::string> node_list(parser, "FILE", "a file with one node id per line", {'l', "node-list"});
        //args::ValueFlag<uint64_t> target_node(parser, "ID", "a single node from which to begin our traversal",{'n', "node"});
        args::ValueFlag<uint64_t> _context_size(parser, "N","the number of steps away from our initial subgraph that we should collect",{'c', "context"});

        /// Range selection
        args::ValueFlag<std::string> _path_range(parser, "STRING",
                                                 "find the node(s) in the specified path range(s) TARGET=path[:pos1[-pos2]]",
                                                 {'p', "path"});

        args::ValueFlag<uint64_t> threads(parser, "N", "number of threads to use", {'t', "threads"});

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

        if (!og_in_file) {
            std::cerr
                    << "[odgi::extract] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
                    << std::endl;
            return 1;
        }

        if (!og_out_file) {
            std::cerr
                    << "[odgi::extract] error: please specify an output file to where to store the extracted graph via -o=[FILE], --out=[FILE]."
                    << std::endl;
            return 1;
        }

        uint64_t context_size = _context_size ? args::get(_context_size) : 0;

        std::vector<string> targets_str;
        if (_path_range) {
            targets_str.push_back(args::get(_path_range));
        }

        graph_t graph;
        assert(argc > 0);
        std::string infile = args::get(og_in_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
                ifstream f(infile.c_str());
                graph.deserialize(f);
                f.close();
            }
        }

        std::vector<Region> targets;

        for (auto &target : targets_str) {
            Region region;
            parse_region(target, region);

            if (!graph.has_path(region.seq)) {
                std::cerr
                        << "[odgi::extract] error: path " << region.seq << " not found in the input graph."
                        << std::endl;
                return 1;
            }
            //std::cerr << region.seq << ": " << region.start << "-" << region.end << std::endl;

            targets.push_back(region);
        }


        if (!targets.empty()) {
            auto prep_graph = [](graph_t& source, graph_t& subgraph, uint64_t context_size) {
                if (context_size > 0) {
                    //if (use_length) {
                    //    algorithms::expand_subgraph_by_length(source, subgraph, context_size);
                    //} else {
                        algorithms::expand_subgraph_by_steps(source, subgraph, context_size, false);
                    //}
                } else {
                    algorithms::add_connecting_edges_to_subgraph(source, subgraph);
                }
                algorithms_add_subpaths_to_subgraph(source, subgraph, false);
                //TODO TODO graph.remove_orphan_edges();
                // Order the mappings by rank. TODO: how do we handle breaks between
                // different sections of a path with a single name?
                //TODO TODO graph.paths.sort_by_mapping_rank();
            };

            for (auto &target : targets) {
                graph_t extract;

                // no coordinates given, we do whole thing (0,-1)
                if (target.start < 0 && target.end < 0) {
                    target.start = 0;
                }

                path_handle_t path_handle = graph.get_path_handle(target.seq);
                algorithms_extract_path_range(graph, path_handle, target.start, target.end, extract);

                prep_graph(graph, extract, context_size);

                std::string outfile = args::get(og_out_file);
                if (!outfile.empty()) {
                    if (outfile == "-") {
                        extract.serialize(std::cout);
                    } else {
                        ofstream f(outfile.c_str());
                        extract.serialize(f);
                        f.close();
                    }
                }
            }

            return 0;
        }

        //todo


        /*
     std::string nucleotide_range = args::get(_path_range);
        if (nucleotide_range.empty()) {
            std::cerr
                    << "[odgi::extract] error: please specify a path nucleotide range: STRING=PATH:start-end."
                    << std::endl;
            return 1;
        }

        size_t foundFirstColon = nucleotide_range.find(':');
        if (foundFirstColon == string::npos) {
            std::cerr
                    << "[odgi::extract] error: please specify a path name."
                    << std::endl;
            return 1;
        }


         std::string path_name = nucleotide_range.substr(0, foundFirstColon);
        nucleotide_range = nucleotide_range.substr(foundFirstColon + 1);

        std::regex regex("-");
        std::vector<std::string> splitted(
                std::sregex_token_iterator(nucleotide_range.begin(), nucleotide_range.end(), regex, -1),
                std::sregex_token_iterator()
        );
        if (splitted.size() != 2) {
            std::cerr
                    << "[odgi::extract] error: please specify a valid nucleotide range: STRING=PATH:start-end."
                    << std::endl;
            return 1;
        }

        if ((splitted[0] != "*" && !utils::is_number(splitted[0])) ||
            (splitted[1] != "*" && !utils::is_number(splitted[1]))) {
            std::cerr
                    << "[odgi::extract] error: please specify valid numbers for the nucleotide range."
                    << std::endl;
            return 1;
        }
        if (!graph.has_path(path_name)) {
            std::cerr
                    << "[odgi::extract] error: please specify a valid path name."
                    << std::endl;
            return 1;
        }

        uint64_t pangenomic_start_pos;
        uint64_t pangenomic_end_pos;

        if (splitted[0] == "*") {
            pangenomic_start_pos = 0;
        } else {
            pangenomic_start_pos = stod(splitted[0]);
        }

        if (splitted[1] == "*") {
            pangenomic_end_pos = std::numeric_limits<uint64_t>::max();
        } else {
            pangenomic_end_pos = stod(splitted[1]);
        }

        if (pangenomic_start_pos > pangenomic_end_pos) {
            std::cerr
                    << "[odgi::extract] error: please specify a start position less than or equal to the end position."
                    << std::endl;
            return 1;
        }

        if (args::get(threads)) {
            omp_set_num_threads(args::get(threads));
        }

        graph_t extract;

        uint64_t new_pangenomic_start_pos = std::numeric_limits<uint64_t>::max();
        uint64_t new_pangenomic_end_pos = 0;

        path_handle_t path_handle = graph.get_path_handle(path_name);

        uint64_t nt_position_in_path = 0;

        graph.for_each_step_in_path(path_handle, [&](const step_handle_t &occ) {
            handle_t h = graph.get_handle_of_step(occ);
            uint64_t h_len = graph.get_length(h);

            //Todo dumb implementation: improve with the math later
            bool take_handle = false;
            do {
                if (nt_position_in_path >= pangenomic_start_pos && nt_position_in_path <= pangenomic_end_pos) {
                    take_handle = true;
                }

                ++nt_position_in_path;
            } while (--h_len != 0);

            if (take_handle) {
                nid_t node_id = graph.get_id(h);
                std::cerr << "node_id: " << node_id << std::endl;

                if (!extract.has_node(node_id)) {
                    if (graph.get_is_reverse(h)) {
                        h = graph.flip(h);
                    }

                    extract.create_handle(graph.get_sequence(h), node_id);
                }
            }
        });

        // expand the context and get path information
        // if forward_only true, then we only go forward.
//            if (context > 0) {
//                algorithms::expand_subgraph_by_steps(*graph, *vg_subgraph, context, forward_only);
//            }
//            if (length > 0) {
//                algorithms::expand_subgraph_by_length(*graph, *vg_subgraph, context, forward_only);
//            } else if (context == 0 && length == 0) {
        algorithms::add_connecting_edges_to_subgraph(graph, extract);
//            }

        algorithms::add_subpaths_to_subgraph(graph, extract);
        */

        // collect the new graph
        /*
        if (args::get(target_node)) {
            uint64_t node_id = args::get(target_node);
            handle_t handle = graph.get_handle(node_id);
            extract.create_handle(graph.get_sequence(handle), node_id);
        }
         */
        graph_t extract;
        if (!args::get(node_list).empty()) {
            ifstream nodes(args::get(node_list));
            std::string next;
            while (std::getline(nodes, next)) {
                uint64_t id = std::stol(next);
                if (graph.has_node(id)) {
                    handle_t h = graph.get_handle(id);
                    extract.create_handle(graph.get_sequence(h), id);
                } else {
                    std::cerr << "[odgi::extract] warning, cannot find node " << id << std::endl;
                }
            }
        }
        /*
if (args::get(context_size)) {
   uint64_t context = args::get(context_size);
   for (uint64_t i = 0; i < context; ++i) {
       // get the edges and connected nodes from the graph to fill out th subgraph
       // for each extract node
       // get its component in the old graph
       // and add it to the extract
       std::vector<handle_t> curr_handles;
       extract.for_each_handle([&](const handle_t& h) {
               curr_handles.push_back(h);
           });
       for (auto& h : curr_handles) {
           handle_t old_h = graph.get_handle(extract.get_id(h));
           graph.follow_edges(old_h, false, [&](const handle_t& c) {
                   handle_t x = extract.create_handle(graph.get_sequence(c), graph.get_id(c));
                   extract.create_edge(h, x);
               });
           graph.follow_edges(old_h, true, [&](const handle_t& c) {
                   handle_t x = extract.create_handle(graph.get_sequence(c), graph.get_id(c));
                   extract.create_edge(x, h);
               });
       }
   }
}
*/

//        std::string outfile = args::get(og_out_file);
//        if (!outfile.empty()) {
//            if (outfile == "-") {
//                extract.serialize(std::cout);
//            } else {
//                ofstream f(outfile.c_str());
//                extract.serialize(f);
//                f.close();
//            }
//        }
        return 0;
    }

    static Subcommand odgi_extract("extract", "extract parts of the graph using paths, positions, and nodes",
                                   PIPELINE, 3, main_extract);


}
