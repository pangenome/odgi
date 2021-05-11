#include "subcommand.hpp"
#include "odgi.hpp"
#include "position.hpp"
#include "args.hxx"
#include "split.hpp"
#include "algorithms/bfs.hpp"
#include <omp.h>

#include "src/algorithms/subgraph/extract.hpp"

namespace odgi {

    using namespace odgi::subcommand;

    int main_depth(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi depth";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser("find the depth of graph as defined by query criteria");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> og_file(parser, "FILE", "compute path depths in this graph", {'i', "input"});

        args::ValueFlag<std::string> _subset_paths(parser, "FILE",
                                                  "compute the depth considering only the paths specified in the FILE; "
                                                  "the file must contain one path name per line and a subset of all paths can be specified.",
                                                  {'s', "subset-paths"});

        args::ValueFlag<std::string> path_name(parser, "PATH_NAME", "compute the depth of the given path in the graph",
                                               {'r', "path"});
        args::ValueFlag<std::string> path_file(parser, "FILE", "compute depth for the paths listed in FILE",
                                               {'R', "paths"});
        args::ValueFlag<std::string> graph_pos(parser, "[node_id][,offset[,(+|-)]*]*",
                                               "compute the depth at the given node, e.g. 7 or 3,4 or 42,10,+ or 302,0,-",
                                               {'g', "graph-pos"});
        args::ValueFlag<std::string> graph_pos_file(parser, "FILE", "a file with one graph position per line",
                                                    {'G', "graph-pos-file"});
        args::ValueFlag<std::string> path_pos(parser, "[path_name][,offset[,(+|-)]*]*",
                                              "return depth at the given path position e.g. chrQ or chr3,42 or chr8,1337,+ or chrZ,3929,-",
                                              {'p', "path-pos"});
        args::ValueFlag<std::string> path_pos_file(parser, "FILE", "a file with one path position per line",
                                                   {'F', "path-pos-file"});
        args::ValueFlag<std::string> bed_input(parser, "FILE", "a BED file of ranges in paths in the graph",
                                               {'b', "bed-input"});

        args::Flag graph_depth(parser, "graph-depth",
                               "compute the depth on each node in the graph",
                               {'d', "graph-depth"});

        args::Flag path_depth(parser, "path-depth",
                              "compute the depth on each base in each path in the graph",
                              {'P', "path-depth"});

        args::Flag summarize_depth(parser, "summarize-graph-depth",
                                   "provide a summary of the depth distribution in the graph",
                                   {'S', "summarize"});

        args::ValueFlag<std::string> _windows_in(parser, "LEN:MIN:MAX",
                                                "write a BED file of path intervals where the depth is between MIN and MAX, "
                                                "merging regions not separated by more than LEN bp",
                                                {'w', "windows-in"});
        args::ValueFlag<std::string> _windows_out(parser, "LEN:MIN:MAX",
                                                 "write a BED file of path intervals where the depth is outside of MIN and MAX, "
                                                 "merging regions not separated by more than LEN bp",
                                                 {'W', "windows-out"});

        args::ValueFlag<uint64_t> _num_threads(parser, "N", "number of threads to use", {'t', "threads"});

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
            std::cerr << "[odgi::depth] error: please specify a target graph via -i=[FILE], --idx=[FILE]." << std::endl;
            return 1;
        }

        if (_windows_in && _windows_out) {
            std::cerr << "[odgi::depth] error: please specify -w/--windows-in or -W/--windows-out, not both." << std::endl;
            return 1;
        }

        uint64_t windows_in_len = 0, windows_in_min = 0, windows_in_max = 0;
        if (_windows_in) {
            if (!algorithms::check_and_get_windows_in_out_parameter(args::get(_windows_in), windows_in_len, windows_in_min, windows_in_max)) {
                std::cerr << "[odgi::depth] error: please specify a valid string (LEN:MIN:MAX) for the -w/--windows-in option." << std::endl;
                return 1;
            }
        }

        uint64_t windows_out_len = 0, windows_out_min = 0, windows_out_max = 0;
        if (_windows_out) {
            if (!algorithms::check_and_get_windows_in_out_parameter(args::get(_windows_out), windows_out_len, windows_out_min, windows_out_max)) {
                std::cerr << "[odgi::depth] error: please specify a valid string (LEN:MIN:MAX) for the -W/--windows-out option." << std::endl;
                return 1;
            }
        }

        odgi::graph_t graph;
        assert(argc > 0);
        if (!args::get(og_file).empty()) {
            const std::string infile = args::get(og_file);
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
                ifstream f(infile.c_str());
                graph.deserialize(f);
                f.close();
            }
        }

        const uint64_t num_threads = args::get(_num_threads) ? args::get(_num_threads) : 1;
        omp_set_num_threads((int) num_threads);

        std::vector<bool> paths_to_consider;
        if (_subset_paths) {
            paths_to_consider.resize(graph.get_path_count() + 1, false);

            std::ifstream refs(args::get(_subset_paths).c_str());
            std::string line;
            while (std::getline(refs, line)) {
                if (!line.empty()) {
                    if (!graph.has_path(line)) {
                        std::cerr << "[odgi::depth] error: path " << line << " not found in graph" << std::endl;
                        exit(1);
                    }

                    paths_to_consider[as_integer(graph.get_path_handle(line))] = true;
                }
            }
        }

        // these options are exclusive (probably we should say with a warning)
        std::vector<odgi::pos_t> graph_positions;
        std::vector<odgi::path_pos_t> path_positions;
        std::vector<odgi::path_range_t> path_ranges;

        auto add_graph_pos = [&graph_positions](const odgi::graph_t &graph,
                                                const std::string &buffer) {
            auto vals = split(buffer, ',');
            /*
            if (vals.size() != 3) {
                std::cerr << "[odgi::depth] error: graph position record is incomplete" << std::endl;
                std::cerr << "[odgi::depth] error: got '" << buffer << "'" << std::endl;
                exit(1); // bail
            }
            */
            uint64_t id = std::stoi(vals[0]);
            if (!graph.has_node(id)) {
                std::cerr << "[odgi::depth] error: no node " << id << " in graph" << std::endl;
                exit(1);
            }
            uint64_t offset = 0;
            if (vals.size() >= 2) {
                offset = std::stoi(vals[1]);
                handle_t h = graph.get_handle(id);
                if (graph.get_length(h) < offset) {
                    std::cerr << "[odgi::depth] error: offset of " << offset << " lies beyond the end of node " << id
                              << std::endl;
                    exit(1);
                }
            }
            bool is_rev = false;
            if (vals.size() == 3) {
                is_rev = vals[2] == "-";
            }
            graph_positions.push_back(make_pos_t(id, is_rev, offset));
        };

        auto add_path_pos = [&path_positions](const odgi::graph_t &graph,
                                              const std::string &buffer) {
            if (!buffer.empty()) {
                auto vals = split(buffer, ',');
                /*
                if (vals.size() != 3) {
                    std::cerr << "[odgi::depth] error: path position record is incomplete" << std::endl;
                    std::cerr << "[odgi::depth] error: got '" << buffer << "'" << std::endl;
                    exit(1); // bail
                }
                */
                auto &path_name = vals[0];
                if (!graph.has_path(path_name)) {
                    std::cerr << "[odgi::depth] error: path " << path_name << " not found in graph" << std::endl;
                    exit(1);
                } else {
                    path_positions.push_back({
                        graph.get_path_handle(path_name),
                        (vals.size() > 1 ? (uint64_t) std::stoi(vals[1]) : 0),
                        (vals.size() == 3 ? vals[2] == "-" : false)
                    });
                }
            }
        };

        auto add_bed_range = [&path_ranges](const odgi::graph_t &graph,
                                            const std::string &buffer) {
            if (!buffer.empty() && buffer[0] != '#') {
                auto vals = split(buffer, '\t');
                /*
                if (vals.size() != 3) {
                    std::cerr << "[odgi::depth] error: path position record is incomplete" << std::endl;
                    std::cerr << "[odgi::depth] error: got '" << buffer << "'" << std::endl;
                    exit(1); // bail
                }
                */
                auto &path_name = vals[0];
                if (!graph.has_path(path_name)) {
                    std::cerr << "[odgi::depth] error: path " << path_name << " not found in graph" << std::endl;
                    exit(1);
                } else {
                    uint64_t start = vals.size() > 1 ? (uint64_t) std::stoi(vals[1]) : 0;
                    uint64_t end = 0;
                    if (vals.size() > 2) {
                        end = (uint64_t) std::stoi(vals[2]);
                    } else {
                        // In the BED format, the end is non-inclusive, unlike start
                        graph.for_each_step_in_path(graph.get_path_handle(path_name), [&](const step_handle_t &s) {
                            end += graph.get_length(graph.get_handle_of_step(s));
                        });
                    }

                    if (start > end) {
                        std::cerr << "[odgi::depth] error: wrong input coordinates in row: " << buffer << std::endl;
                        exit(1);
                    }

                    path_ranges.push_back(
                            {
                                    {
                                            graph.get_path_handle(path_name),
                                            start,
                                            false
                                    },
                                    {
                                            graph.get_path_handle(path_name),
                                            end,
                                            false
                                    },
                                    (vals.size() > 3 && vals[3] == "-"),
                                    buffer
                            });
                }
            }
        };

        if (summarize_depth) {
            // we do nothing here, we iterate over the handles in the graph later
        } else if (graph_depth) {
            graph.for_each_handle([&](const handle_t &h) {
                add_graph_pos(graph, std::to_string(graph.get_id(h)));
            });
        } else if (path_depth) {
            std::vector<path_handle_t> paths;
            graph.for_each_path_handle([&paths](const path_handle_t& path) { paths.push_back(path); });
            // for each path handle
#pragma omp parallel for schedule(dynamic, 1)
            for (auto& path : paths) {
                std::stringstream ss;
                ss << graph.get_path_name(path);
                // for each step
                uint64_t pos = 0;
                graph.for_each_step_in_path(
                    path,
                    [&](const step_handle_t& step) {
                        handle_t handle = graph.get_handle_of_step(step);
                        auto depth = graph.get_step_count(handle) - 1; // don't count self
                        auto next_pos = pos + graph.get_length(handle);
                        while (pos++ < next_pos) {
                            ss << " " << depth;
                        }
                    });
#pragma omp critical (cout)
                std::cout << ss.str() << std::endl;
            }
        } else if (graph_pos) {
            // if we're given a graph_pos, we'll convert it into a path pos
            add_graph_pos(graph, args::get(graph_pos));
        } else if (graph_pos_file) {
            std::ifstream gpos(args::get(graph_pos_file));
            std::string buffer;
            while (std::getline(gpos, buffer)) {
                add_graph_pos(graph, buffer);
            }
        } else if (path_pos) {
            // if given a path pos, we convert it into a path pos in our reference set
            add_path_pos(graph, args::get(path_pos));
        } else if (path_pos_file) {
            // if we're given a file of path depths, we'll convert them all
            std::ifstream refs(args::get(path_pos_file));
            std::string buffer;
            while (std::getline(refs, buffer)) {
                add_path_pos(graph, buffer);
            }
        } else if (bed_input) {
            std::ifstream bed_in(args::get(bed_input));
            std::string buffer;
            while (std::getline(bed_in, buffer)) {
                add_bed_range(graph, buffer);
            }
        } else if (path_name) {
            add_bed_range(graph, args::get(path_name));
        } else if (path_file) {
            // for thing in things
            std::ifstream refs(args::get(path_file));
            std::string line;
            while (std::getline(refs, line)) {
                add_bed_range(graph, line);
            }
        } else if (!_windows_in && !_windows_out){
            // using all the paths in the graph
            graph.for_each_path_handle(
                    [&](const path_handle_t &path) { add_bed_range(graph, graph.get_path_name(path)); });
        }

        auto get_graph_pos = [](const odgi::graph_t &graph,
                                const path_pos_t &pos) {
            const auto path_end = graph.path_end(pos.path);
            uint64_t walked = 0;
            for (step_handle_t s = graph.path_begin(pos.path);
                 s != path_end; s = graph.get_next_step(s)) {
                handle_t h = graph.get_handle_of_step(s);
                uint64_t node_length = graph.get_length(h);
                if (walked + node_length > pos.offset) {
                    return make_pos_t(graph.get_id(h), graph.get_is_reverse(h), pos.offset - walked);
                }
                walked += node_length;
            }

#pragma omp critical (cout)
            std::cerr << "[odgi::depth] warning: position " << graph.get_path_name(pos.path) << ":" << pos.offset
                      << " outside of path" << std::endl;
            return make_pos_t(0, false, 0);
        };

        auto get_offset_in_path = [](const odgi::graph_t &graph,
                                     const path_handle_t &path, const step_handle_t &target) {
            const auto path_end = graph.path_end(path);
            uint64_t walked = 0;
            step_handle_t s = graph.path_begin(path);
            for (; s != target; s = graph.get_next_step(s)) {
                handle_t h = graph.get_handle_of_step(s);
                walked += graph.get_length(h);
            }
            assert(s != path_end);
            return walked;
        };

        auto get_graph_node_depth = [](const odgi::graph_t &graph, const nid_t node_id,
                                          const std::vector<bool>& paths_to_consider) {
            const bool subset_paths = !paths_to_consider.empty();

            uint64_t node_depth = 0;
            std::set<uint64_t> unique_paths;

            const handle_t h = graph.get_handle(node_id);

            graph.for_each_step_on_handle(h, [&](const step_handle_t &occ) {
                if (!subset_paths || paths_to_consider[as_integer(graph.get_path_handle_of_step(occ))]) {
                    ++node_depth;
                    unique_paths.insert(as_integer(graph.get_path(occ)));
                }
            });

            return make_pair(node_depth, unique_paths.size());
        };

        if (_windows_in || _windows_out) {
            std::vector<path_handle_t> paths;
            if (_subset_paths) {
                graph.for_each_path_handle([&](const path_handle_t path) {
                    if (paths_to_consider[as_integer(path)]){
                        paths.push_back(path);
                    }
                });
            } else {
                paths.reserve(graph.get_path_count());
                graph.for_each_path_handle([&](const path_handle_t path) {
                    paths.push_back(path);
                });
            }

            // precompute depths for all handles in parallel
            std::vector<uint64_t> depths(graph.get_node_count() + 1);
            graph.for_each_handle(
                [&](const handle_t& h) {
                    auto id = graph.get_id(h);
                    if (id >= depths.size()) {
                        // require optimized graph to use vector rather than a hash table
                        std::cerr << "[odgi::depth] error: graph is not optimized, apply odgi sort -O" << std::endl;
                        assert(false);
                    }
                    depths[id] = get_graph_node_depth(graph, id, paths_to_consider).first;
                }, true);

            auto in_bounds =
                [&](const handle_t &handle) {
                    uint64_t depth = depths[graph.get_id(handle)];
                    return _windows_in ? (depth >= windows_in_min && depth <= windows_in_max) : (depth < windows_out_min || depth > windows_out_max);
                };

            std::cout << "#path\tstart\tend" << std::endl;

            algorithms::windows_in_out(graph, paths, in_bounds, _windows_in ? windows_in_len : windows_out_len,
                           [&](const std::vector<path_range_t>& path_ranges) {
#pragma omp critical (cout)
                               for (auto path_range : path_ranges) {
                                   std::cout << graph.get_path_name(path_range.begin.path) << "\t"
                                             << path_range.begin.offset << "\t"
                                             << path_range.end.offset << std::endl;
                               }
                           }, num_threads);
        }

        if (summarize_depth) {
            std::cout << "#node.count\tgraph.length\tstep.count\tpath.length\tmean.node.depth\tmean.graph.depth" << std::endl;
            std::atomic<uint64_t> step_count; step_count.store(0);
            std::atomic<uint64_t> node_count; node_count.store(0);
            std::atomic<uint64_t> path_length; path_length.store(0);
            std::atomic<uint64_t> graph_length; graph_length.store(0);
            graph.for_each_handle(
                [&](const handle_t& h) {
                    const nid_t node_id = graph.get_id(h);
                    const auto d = get_graph_node_depth(graph, node_id, paths_to_consider);
                    step_count += d.first;
                    ++node_count;
                    const auto l = graph.get_length(graph.get_handle(node_id));
                    graph_length += l;
                    path_length += l * d.first;
                }, true);
            std::cout << node_count << "\t"
                      << graph_length << "\t"
                      << step_count << "\t"
                      << path_length << "\t"
                      << (double)step_count / (double)node_count << "\t"
                      << (double)path_length / (double)graph_length << std::endl;
        }

        if (!graph_positions.empty()) {
            std::cout << "#node.id\tdepth\tdepth.uniq" << std::endl;
#pragma omp parallel for schedule(dynamic, 1)
            for (auto &pos : graph_positions) {
                const nid_t node_id = id(pos);
                const auto depth = get_graph_node_depth(graph, node_id, paths_to_consider);

#pragma omp critical (cout)
                std::cout << node_id << "\t"
                          << depth.first << "\t"
                          << depth.second << std::endl;
            }
        }

        if (!path_positions.empty()) {
            std::cout << "#path.position\tdepth\tdepth.uniq" << std::endl;
#pragma omp parallel for schedule(dynamic, 1)
            for (auto &path_pos : path_positions) {
                const pos_t pos = get_graph_pos(graph, path_pos);

                const nid_t node_id = id(pos);
                const auto depth = get_graph_node_depth(graph, node_id, paths_to_consider);

#pragma omp critical (cout)
                std::cout << (graph.get_path_name(path_pos.path)) << "," << path_pos.offset << ","
                          << (path_pos.is_rev ? "-" : "+") << "\t"
                          << depth.first << "\t" << depth.second << std::endl;
            }
        }

        if (!path_ranges.empty()) {
            const bool subset_paths = !paths_to_consider.empty();

            std::cout << "#path\tstart\tend\tmean.depth" << std::endl;
#pragma omp parallel for schedule(dynamic, 1)
            for (auto &path_range : path_ranges) {
                const uint64_t start = path_range.begin.offset;
                const uint64_t end = path_range.end.offset;
                const path_handle_t path_handle = path_range.begin.path;

                uint64_t depth = 0;

                uint64_t walked = 0;
                const auto path_end = graph.path_end(path_handle);
                for (step_handle_t cur_step = graph.path_begin(path_handle);
                     cur_step != path_end && walked < end; cur_step = graph.get_next_step(cur_step)) {
                    const handle_t cur_handle = graph.get_handle_of_step(cur_step);
                    const uint64_t cur_length = graph.get_length(cur_handle);
                    walked += cur_length;
                    if (walked >= start) {
                        uint64_t d = 0;
                        if (subset_paths) {
                            graph.for_each_step_on_handle(cur_handle, [&](const step_handle_t &step) {
                                if (paths_to_consider[as_integer(graph.get_path_handle_of_step(step))]) {
                                    ++d;
                                }
                            });

                            if (paths_to_consider[as_integer(path_handle)]) {
                                --d;
                            }
                        } else {
                            // the depth is steps_on_handle - 1
                            d = (graph.get_step_count(cur_handle) - 1);
                        };

                        depth += d * (cur_length
                                           - (walked - cur_length < start ? cur_length - (walked - start) : 0)
                                           - (walked > end ? end - walked : 0));
                    }
                }

#pragma omp critical (cout)
                std::cout << (graph.get_path_name(path_handle)) << "\t" << start << "\t" << end << "\t"
                          << (double) depth / (double) (end - start) << std::endl;
            }
        }

        return 0;
    }

    static Subcommand odgi_depth("depth", "find the depth of graph as defined by query criteria",
                                 PIPELINE, 3, main_depth);

}
