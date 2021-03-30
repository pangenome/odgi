#include "subcommand.hpp"
#include "odgi.hpp"
#include "position.hpp"
#include "args.hxx"
#include "split.hpp"
#include "algorithms/bfs.hpp"
#include <omp.h>

namespace odgi {

    using namespace odgi::subcommand;

    int main_overlap(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi overlap";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser("find the paths touched by the input paths");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> og_file(parser, "FILE", "perform the search in this graph", {'i', "input"});

        args::ValueFlag<std::string> subset_paths(parser, "FILE",
                                                  "perform the search considering only the paths specified in the FILE; "
                                                  "the file must contain one path name per line and a subset of all paths can be specified.",
                                                  {'s', "subset-paths"});

        args::ValueFlag<std::string> path_name(parser, "PATH_NAME", "perform the search of the given path in the graph",
                                               {'r', "path"});
        args::ValueFlag<std::string> path_file(parser, "FILE", "perform the search for the paths listed in FILE",
                                               {'R', "paths"});

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

        if (!og_file || args::get(og_file).empty()) {
            std::cerr << "[odgi::overlap] error: please specify an input graph via -i=[FILE], --idx=[FILE]." << std::endl;
            return 1;
        }

        if ((!path_name || args::get(path_name).empty()) && (!path_file || args::get(path_file).empty())) {
            std::cerr << "[odgi::overlap] error: please specify an input path (-r/--path) or a list of paths (with -R/--paths)." << std::endl;
            return 1;
        }

        odgi::graph_t graph;
        assert(argc > 0);
        std::string infile = args::get(og_file);
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
            ifstream f(infile.c_str());
            graph.deserialize(f);
            f.close();
        }


        uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;
        omp_set_num_threads(num_threads);

        std::vector<path_handle_t> paths_to_consider;
        if (subset_paths && !args::get(subset_paths).empty()) {
            std::ifstream refs(args::get(subset_paths));
            std::string line;
            while (std::getline(refs, line)) {
                if (!line.empty()) {
                    if (!graph.has_path(line)) {
                        std::cerr << "[odgi::overlap] error: path " << line << " not found in graph" << std::endl;
                        exit(1);
                    }

                    paths_to_consider.push_back(graph.get_path_handle(line));
                }
            }
        } else {
            // All paths
            paths_to_consider.reserve(graph.get_path_count());
            graph.for_each_path_handle([&](const path_handle_t path) {
                paths_to_consider.push_back(path);
            });
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
                std::cerr << "[odgi::overlap] error: graph position record is incomplete" << std::endl;
                std::cerr << "[odgi::overlap] error: got '" << buffer << "'" << std::endl;
                exit(1); // bail
            }
            */
            uint64_t id = std::stoi(vals[0]);
            if (!graph.has_node(id)) {
                std::cerr << "[odgi::overlap] error: no node " << id << " in graph" << std::endl;
                exit(1);
            }
            uint64_t offset = 0;
            if (vals.size() >= 2) {
                offset = std::stoi(vals[1]);
                handle_t h = graph.get_handle(id);
                if (graph.get_length(h) < offset) {
                    std::cerr << "[odgi::overlap] error: offset of " << offset << " lies beyond the end of node " << id
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
                    std::cerr << "[odgi::overlap] error: path position record is incomplete" << std::endl;
                    std::cerr << "[odgi::overlap] error: got '" << buffer << "'" << std::endl;
                    exit(1); // bail
                }
                */
                auto &path_name = vals[0];
                if (!graph.has_path(path_name)) {
                    std::cerr << "[odgi::overlap] error: path " << path_name << " not found in graph" << std::endl;
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
                    std::cerr << "[odgi::overlap] error: path position record is incomplete" << std::endl;
                    std::cerr << "[odgi::overlap] error: got '" << buffer << "'" << std::endl;
                    exit(1); // bail
                }
                */
                auto &path_name = vals[0];
                if (!graph.has_path(path_name)) {
                    std::cerr << "[odgi::overlap] error: path " << path_name << " not found in graph" << std::endl;
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
                        std::cerr << "[odgi::overlap] error: wrong input coordinates in row: " << buffer << std::endl;
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

        if (path_name) {
            add_bed_range(graph, args::get(path_name));
        } else {//{ if (path_file) {
            // for thing in things
            std::ifstream refs(args::get(path_file).c_str());
            std::string path_name;
            while (std::getline(refs, path_name)) {
                add_bed_range(graph, path_name);
            }
        }

        auto get_graph_pos = [](const odgi::graph_t &graph,
                                const path_pos_t &pos) {
            auto path_end = graph.path_end(pos.path);
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
            std::cerr << "[odgi::overlap] warning: position " << graph.get_path_name(pos.path) << ":" << pos.offset
                      << " outside of path" << std::endl;
            return make_pos_t(0, false, 0);
        };

        auto get_offset_in_path = [](const odgi::graph_t &graph,
                                     const path_handle_t &path, const step_handle_t &target) {
            auto path_end = graph.path_end(path);
            uint64_t walked = 0;
            step_handle_t s = graph.path_begin(path);
            for (; s != target; s = graph.get_next_step(s)) {
                handle_t h = graph.get_handle_of_step(s);
                walked += graph.get_length(h);
            }
            assert(s != path_end);
            return walked;
        };

        if (!path_ranges.empty()) {
            std::cout << "#path\tpath_touched" << std::endl;

//#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
            for (auto &path_range : path_ranges) {
                // Collect handles crossed by the path in the specified range
                uint64_t start = path_range.begin.offset;
                uint64_t end = path_range.end.offset;
                path_handle_t path_handle = path_range.begin.path;

                std::unordered_set<handle_t> handles;

                uint64_t walked = 0;
                auto path_end = graph.path_end(path_handle);
                for (step_handle_t cur_step = graph.path_begin(path_handle);
                     cur_step != path_end && walked < end; cur_step = graph.get_next_step(cur_step)) {
                    handle_t cur_handle = graph.get_handle_of_step(cur_step);
                    uint64_t cur_length = graph.get_length(cur_handle);
                    walked += cur_length;
                    if (walked >= start) {
                        handles.insert(cur_handle);
                    }
                }

                // Collect paths that cross the collected handles
                std::vector<path_handle_t> touched_path_handles;

                // todo to parallelize?
//#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
                for (path_handle_t p_h : paths_to_consider) {
                    if (p_h != path_handle) {
                        bool stop = false;
                        graph.for_each_step_in_path(p_h, [&](const step_handle_t &step) {
                            handle_t h = graph.get_handle_of_step(step);
                            if (!stop && handles.count(h) > 0) {
                                touched_path_handles.push_back(p_h);
                                stop = true;
                            }
                        });
                    }
                }

//#pragma omp critical (cout)
                for (auto touched_path_handle : touched_path_handles) {
                    std::cout << (graph.get_path_name(path_handle)) << "\t" << start << "\t" << end << "\t"
                              << graph.get_path_name(touched_path_handle) << std::endl;
                }
            }
        }

        return 0;
    }

    static Subcommand odgi_overlap("overlap", "find the paths touched by the input paths",
                                 PIPELINE, 3, main_overlap);

}
