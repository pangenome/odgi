#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include <regex>
#include "utils.hpp"

#include "src/algorithms/subgraph/extract.hpp"

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
        args::Flag _split_subgraphs(parser, "split_subgraphs",
                                    "instead of writing the target subgraphs into a single graph, "
                                    "write one subgraph per given target to a separate file named path:start-end.og  "
                                    "(0-based coordinates)", {'s', "split-subgraphs"});

        args::ValueFlag<uint64_t> _target_node(parser, "ID", "a single node from which to begin our traversal",
                                               {'n', "node"});
        args::ValueFlag<std::string> _node_list(parser, "FILE", "a file with one node id per line", {'l', "node-list"});

        args::ValueFlag<uint64_t> _context_size(parser, "N",
                                                "the number of steps away from our initial subgraph that we should collect",
                                                {'c', "context"});
        args::Flag _use_length(parser, "use_length",
                               "treat the context size as a length in bases (and not as a number of steps)",
                               {'L', "use-length"});

        args::ValueFlag<std::string> _path_range(parser, "STRING",
                                                 "find the node(s) in the specified path range TARGET=path[:pos1[-pos2]] "
                                                 "(0-based coordinates)", {'r', "path-range"});
        args::ValueFlag<std::string> _path_bed_file(parser, "FILE",
                                                    "find the node(s) in the path range(s) specified in the given BED FILE",
                                                    {'b', "bed-file"});
        args::Flag _full_range(parser, "use_length",
                               "collects all nodes in the sorted order of the graph in the min and max position touched by the given path ranges. "
                               "Be careful to use it with very complex graphs",
                               {'E', "full-range"});

        args::ValueFlag<std::string> _path_names_file(parser, "FILE",
                                                      "list of paths to consider in the extraction; the file must "
                                                      "contain one path name per line and a subset of all paths can be specified.",
                                                      {'p', "paths-to-extract"});

        args::ValueFlag<uint64_t> nthreads(parser, "N", "number of threads to use (to embed the subpaths in parallel)",
                                           {'t', "threads"});

        args::Flag _show_progress(parser, "progress",
                                  "print information about the operations and the progress to stderr",
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

        if (!og_in_file) {
            std::cerr
                    << "[odgi::extract] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
                    << std::endl;
            return 1;
        }

        if (_split_subgraphs) {
            if (og_out_file) {
                std::cerr << "[odgi::extract] error: please do not specify an output file (with -o/--out) when "
                             "one subgraph per given target is requested (with -s/--split-subgraphs)." << std::endl;
                return 1;
            }

            if (_target_node) {
                std::cerr << "[odgi::extract] error: please do not specify a single node (with -n/--node) when "
                             "one subgraph per given target is requested (with -s/--split-subgraphs)." << std::endl;
                return 1;
            }

            if (_node_list) {
                std::cerr << "[odgi::extract] error: please do not specify a node list (with -l/--list) when "
                             "one subgraph per given target is requested (with -s/--split-subgraphs)." << std::endl;
                return 1;
            }
        } else {
            if (!og_out_file) {
                std::cerr << "[odgi::extract] error: please specify an output file to where to store the subgraph via "
                             "-o=[FILE], --out=[FILE]." << std::endl;
                return 1;
            }
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

        if (_full_range) {
            nid_t last_node_id = graph.min_node_id();
            graph.for_each_handle([&](const handle_t &h) {
                nid_t node_id = graph.get_id(h);
                if (node_id - last_node_id > 1) {
                    std::cerr
                            << "[odgi::extract] error: the graph is not optimized. "
                               "To extract the full ranges, please run 'odgi sort' using -O/--optimize first"
                            << std::endl;
                    exit(1);
                }
                last_node_id = node_id;
            });
        }


        // Prepare all paths for parallelize the next step (actually, not all paths are always present in the subgraph)
        std::vector<path_handle_t> paths;
        if (args::get(_path_names_file).empty()) {
            paths.reserve(graph.get_path_count());
            graph.for_each_path_handle([&](const path_handle_t path) {
                paths.push_back(path);
            });
        } else {
            std::ifstream path_names_in(args::get(_path_names_file));

            uint64_t num_of_paths_in_file = 0;

            std::vector<bool> path_already_seen;
            path_already_seen.resize(graph.get_path_count(), false);

            std::string line;
            while (std::getline(path_names_in, line)) {
                if (!line.empty()) {
                    if (graph.has_path(line)) {
                        path_handle_t path = graph.get_path_handle(line);
                        uint64_t path_rank = as_integer(path) - 1;
                        if (!path_already_seen[path_rank]) {
                            path_already_seen[path_rank] = true;
                            paths.push_back(path);
                        } else {
                            std::cerr << "[odgi::extract] error: in the path list there are duplicated path names."
                                      << std::endl;
                            exit(1);
                        }
                    }

                    ++num_of_paths_in_file;
                }
            }

            path_names_in.close();

            std::cerr << "[odgi::extract] found " << paths.size() << "/" << num_of_paths_in_file
                      << " paths to consider." << std::endl;

            if (paths.size() == 0) {
                std::cerr << "[odgi::extract] error: no path to consider." << std::endl;
                exit(1);
            }
        }

        std::vector<string> targets_str;
        std::vector<Region> targets;

        // handle targets from BED
        if (_path_bed_file) {
            parse_bed_regions(args::get(_path_bed_file), targets);
        }

        // handle targets from command line
        if (_path_range) {
            targets_str.push_back(args::get(_path_range));
        }

        for (auto &target : targets_str) {
            Region region;
            parse_region(target, region);

            if (!graph.has_path(region.seq)) {
                std::cerr
                        << "[odgi::extract] error: path " << region.seq << " not found in the input graph."
                        << std::endl;
                return 1;
            }

            // no coordinates given, we do whole thing (0,-1)
            if (region.start < 0 && region.end < 0) {
                region.start = 0;
            }

            targets.push_back(region);
        }

        if (_split_subgraphs && targets.empty()) {
            std::cerr << "[odgi::extract] error: please specify at least one target when "
                         "one subgraph per given target is requested (with -s/--split-subgraphs)." << std::endl;
            return 1;
        }

        bool show_progress = args::get(_show_progress);
        uint64_t context_size = _context_size ? args::get(_context_size) : 0;

        uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;
        omp_set_num_threads(num_threads);

        auto prep_graph = [&](graph_t &source, const std::vector<path_handle_t> source_paths, graph_t &subgraph,
                              uint64_t context_size, bool use_length, bool full_range)  {
            if (full_range) {
                // find the start and end node of this and fill things in
                nid_t id_start = std::numeric_limits<nid_t>::max();
                nid_t id_end = 1;
                subgraph.for_each_handle([&](handle_t handle) {
                    nid_t id = subgraph.get_id(handle);
                    id_start = std::min(id_start, id);
                    id_end = std::max(id_end, id);
                });

                algorithms::extract_id_range(source, id_start, id_end, subgraph, show_progress
                                                                                ? "[odgi::extract] collecting all nodes in the path range"
                                                                                : "");
            }

            if (context_size > 0) {
                if (show_progress) {
                    std::cerr << "[odgi::extract] expansion and adding connecting edges" << std::endl;
                }

                if (use_length) {
                    algorithms::expand_subgraph_by_length(source, subgraph, context_size, false,
                                                          show_progress ? "[odgi::extract] adding connecting edges"
                                                                        : "");
                } else {
                    algorithms::expand_subgraph_by_steps(source, subgraph, context_size, false,
                                                         show_progress ? "[odgi::extract] adding connecting edges"
                                                                       : "");
                }
            } else {
                algorithms::add_connecting_edges_to_subgraph(source, subgraph, show_progress
                                                                               ? "[odgi::extract] adding connecting edges"
                                                                               : "");
            }

            algorithms::add_subpaths_to_subgraph(source, paths, subgraph, num_threads,
                                                 show_progress ? "[odgi::extract] adding subpaths" : "");

            // This should not be necessary, if the extraction works correctly
            // graph.remove_orphan_edges();
        };

        auto check_and_create_handle = [&](const graph_t &source, graph_t &subgraph, const nid_t node_id) {
            if (graph.has_node(node_id)) {
                if (!subgraph.has_node(node_id)){
                    subgraph.create_handle(graph.get_sequence(graph.get_handle(node_id)), node_id);
                }
            } else {
                std::cerr << "[odgi::extract] warning, cannot find node " << node_id << std::endl;
            }
        };

        if (_split_subgraphs) {
            for (auto &target : targets) {
                graph_t subgraph;

                path_handle_t path_handle = graph.get_path_handle(target.seq);

                if (show_progress) {
                    std::cerr << "[odgi::extract] extracting path range " << target.seq << ":" << target.start
                              << "-"
                              << target.end << std::endl;
                }
                algorithms::extract_path_range(graph, path_handle, target.start, target.end, subgraph);

                prep_graph(graph, paths, subgraph, context_size, _use_length, _full_range);

                string filename = target.seq + ":" + to_string(target.start) + "-" + to_string(target.end) + ".og";

                if (show_progress) {
                    std::cerr << "[odgi::extract] writing " << filename << std::endl;
                }
                ofstream f(filename);
                subgraph.serialize(f);
                f.close();
            }
        } else {
            graph_t subgraph;

            for (auto &target : targets) {
                path_handle_t path_handle = graph.get_path_handle(target.seq);

                if (show_progress) {
                    std::cerr << "[odgi::extract] extracting path range " << target.seq << ":" << target.start
                              << "-"
                              << target.end << std::endl;
                }
                algorithms::extract_path_range(graph, path_handle, target.start, target.end, subgraph);
            }

            if (args::get(_target_node)) {
                check_and_create_handle(graph, subgraph, args::get(_target_node));
            }

            if (!args::get(_node_list).empty()) {
                ifstream nodes(args::get(_node_list));
                std::string next;
                while (std::getline(nodes, next)) {
                    check_and_create_handle(graph, subgraph, std::stol(next));
                }
            }

            prep_graph(graph, paths, subgraph, context_size, _use_length, _full_range);

            std::string outfile = args::get(og_out_file);
            if (!outfile.empty()) {
                if (outfile == "-") {
                    subgraph.serialize(std::cout);
                } else {
                    ofstream f(outfile.c_str());
                    subgraph.serialize(f);
                    f.close();
                }
            }
        }

        return 0;
    }

    static Subcommand odgi_extract("extract", "extract parts of the graph using paths, positions, and nodes",
                                   PIPELINE, 3, main_extract);

}
