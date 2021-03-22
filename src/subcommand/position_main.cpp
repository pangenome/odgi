#include "subcommand.hpp"
#include "odgi.hpp"
#include "position.hpp"
#include "args.hxx"
#include "split.hpp"
#include "algorithms/bfs.hpp"
#include <omp.h>

namespace odgi {

using namespace odgi::subcommand;

int main_position(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi position";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("position parts of the graph as defined by query criteria");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> og_target_file(parser, "FILE", "describe positions in this graph", {'i', "target"});
    args::ValueFlag<std::string> og_source_file(parser, "FILE", "translate positions from this graph into the target graph using common --lift-paths shared between both graphs [default: use the same source/target graph]", {'x', "source"});
    //args::ValueFlag<std::string> og_out_file(parser, "FILE", "store the graph self index in this file", {'o', "out"});
    args::ValueFlag<std::string> ref_path_name(parser, "PATH_NAME", "translate the given positions into positions relative to this reference path", {'r', "ref-path"});
    args::ValueFlag<std::string> ref_path_file(parser, "FILE", "use the ref-paths in FILE", {'R', "ref-paths"});
    args::ValueFlag<std::string> lift_path_name(parser, "PATH_NAME", "lift positions from --source to --target via coordinates in this path common to both graphs [default: all common paths between --source and --target]", {'l', "lift-path"});
    args::ValueFlag<std::string> lift_path_file(parser, "FILE", "use the lift-paths in FILE", {'L', "lift-paths"});
    args::ValueFlag<std::string> graph_pos(parser, "[node_id][,offset[,(+|-)]*]*", "a graph position, e.g. 42,10,+ or 302,0,-", {'g', "graph-pos"});
    args::ValueFlag<std::string> graph_pos_file(parser, "FILE", "a file with one graph position per line", {'G', "graph-pos-file"});
    args::ValueFlag<std::string> path_pos(parser, "[path_name][,offset[,(+|-)]*]*", "a path position, e.g. chr8,1337,+ or chrZ,3929,-", {'p', "path-pos"});
    args::ValueFlag<std::string> path_pos_file(parser, "FILE", "a file with one path position per line", {'F', "path-pos-file"});
    args::ValueFlag<std::string> bed_input(parser, "FILE", "a BED file of ranges in paths in the graph to lift into the target graph", {'b', "bed-input"});
    args::Flag give_graph_pos(parser, "give-graph-pos", "emit graph positions (node,offset,strand) rather than path positions", {'v', "give-graph-pos"});
    args::ValueFlag<uint64_t> _search_radius(parser, "DISTANCE", "limit coordinate conversion breadth-first search up to DISTANCE bp from each given position [default: 10000]", {'d',"search-radius"});
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
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    if (!og_target_file) {
        std::cerr << "[odgi::position] error: please specify a target graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    odgi::graph_t target_graph;
    assert(argc > 0);
    std::string infile = args::get(og_target_file);
    if (infile.size()) {
        if (infile == "-") {
            target_graph.deserialize(std::cin);
        } else {
            ifstream f(infile.c_str());
            target_graph.deserialize(f);
            f.close();
        }
    }

    bool lifting = false;
    odgi::graph_t source_graph; // will be empty if we don't have different source and target
    if (og_source_file) {
        lifting = true;
        std::string infile = args::get(og_source_file);
        if (infile.size()) {
            if (infile == "-") {
                source_graph.deserialize(std::cin);
            } else {
                ifstream f(infile.c_str());
                source_graph.deserialize(f);
                f.close();
            }
        }
    }

    if (args::get(threads)) {
        omp_set_num_threads(args::get(threads));
    }

    // todo: load many positions from a file
    // todo: convert a BED file
    // to simplify parallelism, collect our positions when doing so

    // collect our reference paths
    std::vector<path_handle_t> ref_paths;
    if (ref_path_name) {
        std::string path_name = args::get(ref_path_name);
        if (!target_graph.has_path(path_name)) {
            std::cerr << "[odgi::position] error: ref path " << path_name << " not found in graph" << std::endl;
            return 1;
        } else {
            ref_paths.push_back(target_graph.get_path_handle(path_name));
        }
    } else if (ref_path_file) {
        // for thing in things
        std::ifstream refs(args::get(ref_path_file).c_str());
        std::string path_name;
        while (std::getline(refs, path_name)) {
            if (!target_graph.has_path(path_name)) {
                std::cerr << "[odgi::position] error: ref path " << path_name << " not found in graph" << std::endl;
                return 1;
            } else {
                ref_paths.push_back(target_graph.get_path_handle(path_name));
            }
        }
    } else {
        // using all the paths in the graph
        target_graph.for_each_path_handle([&](const path_handle_t& path) { ref_paths.push_back(path); });
    }

    // for translating source to target
    std::vector<path_handle_t> lift_paths_source;
    std::vector<path_handle_t> lift_paths_target;
    if ((lift_path_name || lift_path_file) && !lifting) {
        std::cerr << "[odgi::position] error: lifting requires a separate source and target graph, specify --source" << std::endl;
        return 1;
    } else if (lifting) {
        if (lift_path_name) {
            std::string path_name = args::get(lift_path_name);
            if (!target_graph.has_path(path_name)
                || !source_graph.has_path(path_name)) {
                std::cerr << "[odgi::position] error: lift path " << path_name << " not found in both source and target graph" << std::endl;
                return 1;
            } else {
                lift_paths_source.push_back(source_graph.get_path_handle(path_name));
                lift_paths_target.push_back(target_graph.get_path_handle(path_name));
            }
        } else if (lift_path_file) {
            // for thing in things
            std::ifstream refs(args::get(lift_path_file).c_str());
            std::string path_name;
            while (std::getline(refs, path_name)) {
                if (!target_graph.has_path(path_name)
                    || !source_graph.has_path(path_name)) {
                    std::cerr << "[odgi::position] error: lift path " << path_name << " not found in both source and target graph" << std::endl;
                    return 1;
                } else {
                    lift_paths_source.push_back(source_graph.get_path_handle(path_name));
                    lift_paths_target.push_back(target_graph.get_path_handle(path_name));
                }
            }
        } else {
            // using the common set of paths between the two graphs
            // make a set intersection
            std::vector<std::string> lift_names_source, lift_names_target, lift_names_common;
            source_graph.for_each_path_handle([&](const path_handle_t& path) {
                                                  lift_names_source.push_back(source_graph.get_path_name(path)); });
            target_graph.for_each_path_handle([&](const path_handle_t& path) {
                                                  lift_names_target.push_back(target_graph.get_path_name(path)); });
            std::sort(lift_names_source.begin(), lift_names_source.end());
            std::sort(lift_names_target.begin(), lift_names_target.end());
            std::set_intersection(lift_names_source.begin(), lift_names_source.end(),
                                  lift_names_target.begin(), lift_names_target.end(),
                                  std::back_inserter(lift_names_common));
            for (auto& path_name : lift_names_common) {
                lift_paths_source.push_back(source_graph.get_path_handle(path_name));
                lift_paths_target.push_back(target_graph.get_path_handle(path_name));
            }
            //target_graph.
            //lift_paths
        }
        if (lift_paths_source.size() != lift_paths_target.size()) {
            std::cerr << "[odgi::position] error: differing number of lift paths in target and source, suggests error" << std::endl;
            return 1;
        }
        if (lift_paths_source.empty() || lift_paths_target.empty()) {
            std::cerr << "[odgi::position] error: no lift paths common to both target and source, cannot proceed" << std::endl;
            std::cerr << "[odgi::position] select a set of common paths as lift paths or ensure that there are paths in common" << std::endl;
            return 1;
        }
    }

    // these options are exclusive (probably we should say with a warning)
    std::vector<odgi::pos_t> graph_positions;
    std::vector<odgi::path_pos_t> path_positions;
    std::vector<odgi::path_range_t> path_ranges;

    // TODO the parsers here should find the last 2 , delimiters and split on them
    
    auto add_graph_pos =
        [&graph_positions](const odgi::graph_t& graph,
                           const std::string& buffer) {
            auto vals = split(buffer, ',');
            /*
            if (vals.size() != 3) {
                std::cerr << "[odgi::position] error: graph position record is incomplete" << std::endl;
                std::cerr << "[odgi::position] error: got '" << buffer << "'" << std::endl;
                exit(1); // bail
            }
            */
            uint64_t id = std::stoi(vals[0]);
            if (!graph.has_node(id)) {
                std::cerr << "[odgi::position] error: no node " << id << " in graph" << std::endl;
                exit(1);
            }
            uint64_t offset = 0;
            if (vals.size() >= 2) {
                offset = std::stoi(vals[1]);
                handle_t h = graph.get_handle(id);
                if (graph.get_length(h) < offset) {
                    std::cerr << "[odgi::position] error: offset of " << offset << " lies beyond the end of node " << id << std::endl;
                    exit(1);
                }
            }
            bool is_rev = false;
            if (vals.size() == 3) {
                is_rev = vals[2] == "-";
            }
            graph_positions.push_back(make_pos_t(id, is_rev, offset));
        };

    auto add_path_pos =
        [&path_positions](const odgi::graph_t& graph,
                          const std::string& buffer) {
            auto vals = split(buffer, ',');
            /*
            if (vals.size() != 3) {
                std::cerr << "[odgi::position] error: path position record is incomplete" << std::endl;
                std::cerr << "[odgi::position] error: got '" << buffer << "'" << std::endl;
                exit(1); // bail
            }
            */
            auto& path_name = vals[0];
            if (!graph.has_path(path_name)) {
                std::cerr << "[odgi::position] error: ref path " << path_name << " not found in graph" << std::endl;
                exit(1);
            } else {
                path_positions.push_back({
                        graph.get_path_handle(path_name),
                        (vals.size() > 1 ? (uint64_t)std::stoi(vals[1]) : 0),
                        (vals.size() == 3 ? vals[2] == "-" : false)
                    });
            }
        };

    auto add_bed_range =
        [&path_ranges](const odgi::graph_t& graph,
                       const std::string& buffer) {
            if (!buffer.empty() && buffer[0] != '#') {
                auto vals = split(buffer, '\t');
                /*
                if (vals.size() != 3) {
                    std::cerr << "[odgi::position] error: path position record is incomplete" << std::endl;
                    std::cerr << "[odgi::position] error: got '" << buffer << "'" << std::endl;
                    exit(1); // bail
                }
                */
                auto& path_name = vals[0];
                if (!graph.has_path(path_name)) {
                    std::cerr << "[odgi::position] error: ref path " << path_name << " not found in graph" << std::endl;
                    exit(1);
                } else {
                    uint64_t start = vals.size() > 1 ? (uint64_t) std::stoi(vals[1]) : 0;
                    uint64_t end = 0;
                    if (vals.size() > 2) {
                        end = (uint64_t) std::stoi(vals[2]);
                    } else {
                        graph.for_each_step_in_path(graph.get_path_handle(path_name), [&](const step_handle_t &s) {
                            end += graph.get_length(graph.get_handle_of_step(s));
                        });
                        --end; // BED-style to 0-based inclusive coordinates
                    }

                    if (start > end) {
                        std::cerr << "[odgi::position] error: wrong input coordinates in row: " << buffer << std::endl;
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
    
    if (graph_pos) {
        // if we're given a graph_pos, we'll convert it into a path pos
        if (lifting) {
            add_graph_pos(source_graph, args::get(graph_pos));
        } else {
            add_graph_pos(target_graph, args::get(graph_pos));
        }
    } else if (graph_pos_file) {
        std::ifstream gpos(args::get(graph_pos_file).c_str());
        std::string buffer;
        while (std::getline(gpos, buffer)) {
            add_graph_pos(source_graph, buffer);
        }
    } else if (path_pos) {
        // if given a path pos, we convert it into a path pos in our reference set
        if (lifting) {
            add_path_pos(source_graph, args::get(path_pos));
        } else {
            add_path_pos(target_graph, args::get(path_pos));
        }
    } else if (path_pos_file) {
        // if we're given a file of path positions, we'll convert them all
        std::ifstream refs(args::get(path_pos_file).c_str());
        std::string buffer;
        while (std::getline(refs, buffer)) {
            if (lifting) {
                add_path_pos(source_graph, buffer);
            } else {
                add_path_pos(target_graph, buffer);
            }
        }
    } else if (bed_input) {
        std::ifstream bed_in(args::get(bed_input).c_str());
        std::string buffer;
        while (std::getline(bed_in, buffer)) {
            if (lifting) {
                add_bed_range(source_graph, buffer);
            } else {
                add_bed_range(target_graph, buffer);
            }
        }
    }
    // todo: bed files

    uint64_t search_radius = _search_radius ? args::get(_search_radius) : 10000;

    // make an hash set of our ref path ids for quicker lookup
    hash_set<uint64_t> ref_path_set;
    for (auto& path : ref_paths) {
        ref_path_set.insert(as_integer(path));
    }
    hash_set<uint64_t> lift_path_set_source;
    hash_set<uint64_t> lift_path_set_target;
    for (auto& path : lift_paths_source) {
        lift_path_set_source.insert(as_integer(path));
    }
    for (auto& path : lift_paths_target) {
        lift_path_set_target.insert(as_integer(path));
    }

    auto get_graph_pos =
        [](const odgi::graph_t& graph,
           const path_pos_t& pos) {
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
            std::cerr << "[odgi::position] warning: position " << graph.get_path_name(pos.path) << ":" << pos.offset << " outside of path" << std::endl;
            return make_pos_t(0, false, 0);
        };

    auto get_offset_in_path =
        [](const odgi::graph_t& graph,
           const path_handle_t& path, const step_handle_t& target) {
            auto path_end = graph.path_end(path);
            uint64_t walked = 0;
            step_handle_t s = graph.path_begin(path);
            for ( ;  s != target; s = graph.get_next_step(s)) {
                handle_t h = graph.get_handle_of_step(s);
                walked += graph.get_length(h);
            }
            assert(s != path_end);
            return walked;
        };

    // TODO this part needs to include adjustments for in-node offsets vs. where we find the ref path
    // TODO should we always look "backwards" when seeking the ref pos?

    struct lift_result_t {
        int64_t path_offset = 0;
        step_handle_t ref_hit;
        uint64_t walked_to_hit_ref = 0;
        bool is_rev_vs_ref = false;
        bool used_bidirectional = false;
    };

    auto get_position =
        [&search_radius,&get_offset_in_path](const odgi::graph_t& graph,
                                             const hash_set<uint64_t>& path_set,
                                             const pos_t& pos, lift_result_t& lift) {
            // unpacking our args
            int64_t& path_offset = lift.path_offset;
            step_handle_t& ref_hit = lift.ref_hit;
            uint64_t& walked_to_hit_ref = lift.walked_to_hit_ref;
            bool& rev_vs_ref = lift.is_rev_vs_ref;
            bool& used_bidirectional = lift.used_bidirectional;
            handle_t start_handle = graph.get_handle(id(pos), is_rev(pos));
            bool found_hit = false;
            uint64_t adj_last_node = 0;
            hash_set<uint64_t> seen;
            for (auto try_bidirectional : { false, true }) {
                if (try_bidirectional) used_bidirectional = true;
                odgi::algorithms::bfs(
                    graph,
                    [&](const handle_t& h, const uint64_t& r, const uint64_t& l, const uint64_t& d) {
                        seen.insert(as_integer(h));
                        bool got_hit = false;
                        step_handle_t hit;
                        graph.for_each_step_on_handle(
                            h, [&](const step_handle_t& s) {
                                   auto p = graph.get_path_handle_of_step(s);
                                   if (!got_hit && path_set.count(as_integer(p))) {
                                       //std::cerr << "thought I got a hit" << std::endl;
                                       got_hit = true;
                                       hit = s;
                                       walked_to_hit_ref += l; // how far we came to get to this node
                                       rev_vs_ref = graph.get_is_reverse(graph.get_handle_of_step(s)) == graph.get_is_reverse(h);
                                       if (d == 0) { // if we're on the start node
                                           if (rev_vs_ref) {
                                               // and if the path orientation is the same as our traversal orientation
                                               // then we need to add the remaining distance from our original offset to the end of node
                                               // to the final path position offset
                                               adj_last_node = graph.get_length(h) - offset(pos);
                                           } else {
                                               // otherwise if the original path is in the same orientation
                                               // then we add the original forward offset to the ref path offset
                                               adj_last_node = offset(pos);
                                           }
                                       } else { // if we're not on the first node
                                           if (rev_vs_ref) {
                                               // and we come onto the result in the same orientation
                                               // it means the ref pos is at the node end
                                               adj_last_node = 0; // so we have no adjustment
                                           } else {
                                               // otherwise, it means the original path is in the same orientation
                                               // then we need to adjust by the length of this stepb
                                               // because we enter at node end, but we'll get the graph position for the step
                                               // at the node beginning
                                               adj_last_node = graph.get_length(h);
                                           }
                                       }
                                   }
                               });
                        if (got_hit) {
                            ref_hit = hit;
                            found_hit = true;
                        }
                    },
                    [&seen](const handle_t& h) { return seen.count(as_integer(h)); },
                    [](const handle_t& l, const handle_t& h) { return false; },
                    [&found_hit](void) { return found_hit; },
                    { graph.flip(start_handle) },
                    { },
                    try_bidirectional,
                    0,
                    search_radius);
                if (found_hit) break; // if we got a hit, don't go bidirectional
            }
            if (found_hit) {
                path_handle_t p = graph.get_path_handle_of_step(ref_hit);
                // TODO ORIENTATION
                path_offset = get_offset_in_path(graph, p, ref_hit) + adj_last_node;
                return true;
            } else {
                path_offset = -1;
                return false;
            }
        };

    // for each position that we want to look up
#pragma omp parallel for schedule(dynamic,1)
    for (auto& _pos : graph_positions) {
        // go to the graph
        // do a little BFS, bounded by our limit
        // now, if we found our hit, print
        // optionally, we will translate from a source graph into a target graph
        lift_result_t source_result;
        //bool ok = false;
        pos_t pos;
        if (lifting) {
            if (get_position(source_graph, lift_path_set_source, _pos, source_result)) {
                pos = get_graph_pos(target_graph,
                                    { target_graph.get_path_handle(
                                            source_graph.get_path_name(
                                                source_graph.get_path_handle_of_step(
                                                    source_result.ref_hit))),
                                      (uint64_t)source_result.path_offset,
                                      source_result.is_rev_vs_ref });
            } else {
                pos = make_pos_t(0,false,0); // couldn't lift
            }
        } else {
            pos = _pos;
        }
        //= (lifting ? ok = get_position(source_graph, lay_path_set
        lift_result_t result;
        if (id(pos) && give_graph_pos) {
            // force graph position in target
#pragma omp critical (cout)
            std::cout << "#source.path.pos\ttarget.graph.pos" << std::endl
                      << id(_pos) << "," << offset(_pos) << "," << (is_rev(_pos) ? "-" : "+") << "\t"
                      << "\t" << id(pos) << "," << offset(pos) << "," << (is_rev(pos) ? "-" : "+") << std::endl;
        } else if (get_position(target_graph, ref_path_set, pos, result)) {
            bool ref_is_rev = false;
            path_handle_t p = target_graph.get_path_handle_of_step(result.ref_hit);
#pragma omp critical (cout)
            std::cout << "#source.path.pos\ttarget.path.pos\tdist.to.ref\tstrand.vs.ref" << std::endl
                      << id(pos) << "," << offset(pos) << "," << (is_rev(pos) ? "-" : "+") << "\t"
                      << target_graph.get_path_name(p) << "," << result.path_offset << "," << (ref_is_rev ? "-" : "+") << "\t"
                      << result.walked_to_hit_ref << "\t" << (result.is_rev_vs_ref ? "-" : "+") << std::endl;
        }
    }

#pragma omp parallel for schedule(dynamic,1)
    for (auto& path_pos : path_positions) {
        // TODO we need a better input format
        pos_t pos;
        // handle the lift into the target graph
        if (lifting) {
            lift_result_t source_result;
            pos_t _pos = get_graph_pos(source_graph, path_pos);
            if (id(_pos) && get_position(source_graph, lift_path_set_source, _pos, source_result)) {
                pos = get_graph_pos(target_graph,
                                    { target_graph.get_path_handle(
                                            source_graph.get_path_name(
                                                source_graph.get_path_handle_of_step(
                                                    source_result.ref_hit))),
                                      (uint64_t)source_result.path_offset,
                                      source_result.is_rev_vs_ref });
            } else {
                pos = make_pos_t(0,false,0); // couldn't lift
            }

        } else {
            //path_pos = _path_pos;
            pos = get_graph_pos(target_graph, path_pos);
        }
        lift_result_t result;
        //std::cerr << "Got graph pos " << id(pos) << std::endl;
        if (id(pos)) {
            if (give_graph_pos) {
#pragma omp critical (cout)
                std::cout << "#source.path.pos\ttarget.graph.pos" << std::endl
                          << (lifting ? source_graph.get_path_name(path_pos.path) : target_graph.get_path_name(path_pos.path))
                          << "," << path_pos.offset << "," << (path_pos.is_rev ? "-" : "+")
                          << "\t" << id(pos) << "," << offset(pos) << "," << (is_rev(pos) ? "-" : "+") << std::endl;
            } else if (get_position(target_graph, ref_path_set, pos, result)) {
                bool ref_is_rev = false;
                path_handle_t p = target_graph.get_path_handle_of_step(result.ref_hit);
#pragma omp critical (cout)
                std::cout << "#source.path.pos\ttarget.path.pos\tdist.to.ref\tstrand.vs.ref" << std::endl
                          << (lifting ? source_graph.get_path_name(path_pos.path) : target_graph.get_path_name(path_pos.path)) << ","
                          << path_pos.offset << "," << (path_pos.is_rev ? "-" : "+") << "\t"
                          << target_graph.get_path_name(p) << "," << result.path_offset << "," << (ref_is_rev ? "-" : "+") << "\t"
                          << result.walked_to_hit_ref << "\t" << (result.is_rev_vs_ref ? "-" : "+") << std::endl;
            }
        }
    }

#pragma omp parallel for schedule(dynamic,1)
    for (auto& path_range : path_ranges) {
        pos_t pos_begin, pos_end;
        // handle the lift into the target graph
        if (lifting) {
            lift_result_t source_begin_result, source_end_result;
            pos_t _pos_begin = get_graph_pos(source_graph, path_range.begin);
            pos_t _pos_end = get_graph_pos(source_graph, path_range.end);
            if (id(_pos_begin) && get_position(source_graph, lift_path_set_source, _pos_begin, source_begin_result)
                && id(_pos_end) && get_position(source_graph, lift_path_set_source, _pos_end, source_end_result)) {
                pos_begin = get_graph_pos(target_graph,
                                          { target_graph.get_path_handle(
                                                  source_graph.get_path_name(
                                                      source_graph.get_path_handle_of_step(
                                                          source_begin_result.ref_hit))),
                                            (uint64_t)source_begin_result.path_offset,
                                            source_begin_result.is_rev_vs_ref });
                pos_end = get_graph_pos(target_graph,
                                        { target_graph.get_path_handle(
                                                source_graph.get_path_name(
                                                    source_graph.get_path_handle_of_step(
                                                        source_end_result.ref_hit))),
                                          (uint64_t)source_end_result.path_offset,
                                          source_end_result.is_rev_vs_ref });
            } else {
                pos_begin = make_pos_t(0,false,0); // couldn't lift
                pos_end = make_pos_t(0,false,0); // couldn't lift
            }
        } else {
            //path_pos = _path_pos;
            pos_begin = get_graph_pos(target_graph, path_range.begin);
            pos_end = get_graph_pos(target_graph, path_range.end);
        }
        if (id(pos_begin) && id(pos_end)) {
            lift_result_t lift_begin;
            lift_result_t lift_end;
            // TODO add a GAF-style path to the record to say where the BED range walks in the graph
            // TODO optionally list out the nodes in this particular range (e.g. those within it in our sort order)
            if (give_graph_pos) {
#pragma omp critical (cout)
                std::cout << path_range.data << "\t"
                          << id(pos_begin) << "," << offset(pos_begin) << "," << (is_rev(pos_begin)?"-":"+") << "\t"
                          << id(pos_end) << "," << offset(pos_end) << "," << (is_rev(pos_end)?"-":"+") << std::endl;
            } else if (get_position(target_graph, ref_path_set, pos_begin, lift_begin)
                       && get_position(target_graph, ref_path_set, pos_end, lift_end)) {
                bool ref_is_rev = false;
                path_handle_t p_begin = target_graph.get_path_handle_of_step(lift_begin.ref_hit);
                path_handle_t p_end = target_graph.get_path_handle_of_step(lift_end.ref_hit);
                // XXX TODO assert these to be equal......
#pragma omp critical (cout)
                std::cout << path_range.data << "\t"
                          << target_graph.get_path_name(p_begin) << ","
                          << lift_begin.path_offset << ","
                          << (lift_begin.is_rev_vs_ref ? "-" : "+") << "\t"
                          << target_graph.get_path_name(p_end) << ","
                          << lift_end.path_offset << ","
                          << (lift_end.is_rev_vs_ref ? "-" : "+") << "\t"
                          << (lift_begin.is_rev_vs_ref ^ path_range.is_rev ? "-" : "+") << std::endl;
                    //<< walked_to_hit_ref << "\t" << (is_rev_vs_ref ? "-" : "+") << std::endl;
            }
        }
    }

    // todo - lift the position into another graph
    // requires an input of target paths in the final graph
    // and optionally the set of paths in common (we can compute this by default) to drive the lift

    return 0;
}

static Subcommand odgi_position("position", "coordinate projections between nodes and paths",
                                PIPELINE, 3, main_position);


}
