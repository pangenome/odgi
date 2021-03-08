#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include <regex>
#include "utils.hpp"
#include "algorithms/explode.hpp"

namespace odgi {

    using namespace odgi::subcommand;

    int main_extract(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi extract";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser("extract parts of the graph as defined by query criteria");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> og_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
        args::ValueFlag<std::string> og_out_file(parser, "FILE", "store the graph self index in this file",
                                                 {'o', "out"});

        /// Range selection
        args::ValueFlag<std::string> _nucleotide_range(parser, "STRING",
                                                       "nucleotide range to extract: STRING=PATH:start-end. `*-end` for `[0,end]`; `start-*` for `[start,pangenome_length]`.",
                                                       {'r', "path-range"});

        //args::ValueFlag<std::string> node_list(parser, "FILE", "a file with one node id per line", {'l', "node-list"});
        //args::ValueFlag<uint64_t> target_node(parser, "ID", "a single node from which to begin our traversal", {'n', "node"});
        //args::ValueFlag<uint64_t> context_size(parser, "N", "the number of steps away from our initial subgraph that we should collect", {'c', "context"});

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

        // todo... for explode
        if (!og_out_file) {
            std::cerr
                    << "[odgi::extract] error: please specify an output file to where to store the extracted graph via -o=[FILE], --out=[FILE]."
                    << std::endl;
            return 1;
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
        if (args::get(threads)) {
            omp_set_num_threads(args::get(threads));
        }

        graph_t extract;

        std::string nucleotide_range = args::get(_nucleotide_range);
        if (!nucleotide_range.empty()) {
            size_t foundFirstColon = nucleotide_range.find(':');
            if (foundFirstColon == string::npos) {
                std::cerr
                        << "[odgi::extract] error: please specify a path name."
                        << std::endl;
                return 1;
            }

            std::string path_name = nucleotide_range.substr(0, foundFirstColon);
            if (!graph.has_path(path_name)) {
                std::cerr
                        << "[odgi::extract] error: please specify a valid path name."
                        << std::endl;
                return 1;
            }

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

            if (pangenomic_start_pos >= pangenomic_end_pos) {
                std::cerr
                        << "[odgi::extract] error: please specify a start position less than the end position."
                        << std::endl;
                return 1;
            }

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
//        if (context > 0) {
//            algorithms::expand_subgraph_by_steps(*graph, *vg_subgraph, context, forward_only);
//        }
//        if (length > 0) {
//            algorithms::expand_subgraph_by_length(*graph, *vg_subgraph, context, forward_only);
//        }
//        else if (context == 0 && length == 0) {
            algorithms::add_connecting_edges_to_subgraph(graph, extract);
//        }
            //algorithms::add_subpaths_to_subgraph(graph, extract, true);
        }



        // collect the new graph
        /*
        if (args::get(target_node)) {
            uint64_t node_id = args::get(target_node);
            handle_t handle = graph.get_handle(node_id);
            extract.create_handle(graph.get_sequence(handle), node_id);
        }
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
        return 0;
    }

    static Subcommand odgi_extract("extract", "extract parts of the graph using paths, positions, and nodes",
                                   PIPELINE, 3, main_extract);


}
