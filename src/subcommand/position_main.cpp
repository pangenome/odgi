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
    args::ValueFlag<std::string> og_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> og_out_file(parser, "FILE", "store the graph self index in this file", {'o', "out"});
    args::ValueFlag<std::string> graph_pos(parser, "[node_id]:[offset]:[+|-]", "a graph position, e.g. 42:10:+ or 302:0:-", {'g', "graph-pos"});
    args::ValueFlag<std::string> graph_pos_file(parser, "FILE", "a file with one graph position per line", {'G', "graph-pos-file"});
    args::ValueFlag<std::string> path_pos(parser, "[path_name]:[offset]:[+|-]", "a path position, e.g. chr8:1337:+ or chrZ:3929:-", {'p', "path-pos"});
    args::ValueFlag<std::string> path_pos_file(parser, "FILE", "a file with one path position per line", {'F', "path-pos-file"});
    args::ValueFlag<std::string> ref_path_name(parser, "PATH_NAME", "translate the given positions into positions relative to this reference path", {'r', "ref-path"});
    args::ValueFlag<std::string> ref_path_file(parser, "FILE", "translate the given positions into the positions in paths named in this file", {'R', "ref-paths"});
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

    if (!og_in_file) {
        std::cerr << "[odgi position] error: please specify an input file via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    graph_t graph;
    assert(argc > 0);
    std::string infile = args::get(og_in_file);
    if (infile.size()) {
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

    // todo: load many positions from a file
    // todo: convert a BED file
    // to simplify parallelism, collect our positions when doing so

    // collect our reference paths
    std::vector<path_handle_t> ref_paths;
    if (ref_path_name) {
        std::string path_name = args::get(ref_path_name);
        if (!graph.has_path(path_name)) {
            std::cerr << "[odgi position] error: ref path " << path_name << " not found in graph" << std::endl;
            return 1;
        } else {
            ref_paths.push_back(graph.get_path_handle(path_name));
        }
    } else if (ref_path_file) {
        // for thing in things
        std::ifstream refs(args::get(ref_path_file).c_str());
        std::string path_name;
        while (std::getline(refs, path_name)) {
            if (!graph.has_path(path_name)) {
                std::cerr << "[odgi position] error: ref path " << path_name << " not found in graph" << std::endl;
                return 1;
            } else {
                ref_paths.push_back(graph.get_path_handle(path_name));
            }
        }
    } else {
        // using all the paths in the graph
        graph.for_each_path_handle([&](const path_handle_t& path) { ref_paths.push_back(path); });
    }

    // these options are exclusive (probably we should say with a warning)
    std::vector<odgi::pos_t> graph_positions;
    std::vector<odgi::path_pos_t> path_positions;

    // TODO the parsers here should find the last 2 , delimiters and split on them
    
    auto add_graph_pos =
        [&](const std::string& buffer) {
            auto vals = split(buffer, ',');
            if (vals.size() != 3) {
                std::cerr << "[odgi position] error: graph position record is incomplete" << std::endl;
                std::cerr << "[odgi position] error: got '" << buffer << "'" << std::endl;
                exit(1); // bail
            }
            uint64_t id = std::stoi(vals[0]);
            uint64_t offset = std::stoi(vals[1]);
            bool is_rev = vals[2] == "-";
            if (!graph.has_node(id)) {
                std::cerr << "[odgi position] error: no node " << id << " in graph" << std::endl;
                exit(1);
            }
            handle_t h = graph.get_handle(id);
            if (graph.get_length(h) < offset) {
                std::cerr << "[odgi position] error: offset of " << offset << " lies beyond the end of node " << id << std::endl;
                exit(1);
            }
            graph_positions.push_back(make_pos_t(id, is_rev, offset));
        };

    auto add_path_pos =
        [&](const std::string& buffer) {
            auto vals = split(buffer, ',');
            if (vals.size() != 3) {
                std::cerr << "[odgi position] error: path position record is incomplete" << std::endl;
                std::cerr << "[odgi position] error: got '" << buffer << "'" << std::endl;
                exit(1); // bail
            }
            auto& path_name = vals[0];
            if (!graph.has_path(path_name)) {
                std::cerr << "[odgi position] error: ref path " << path_name << " not found in graph" << std::endl;
                exit(1);
            } else {
                path_positions.push_back({
                        graph.get_path_handle(path_name),
                        (uint64_t)std::stoi(vals[1]),
                        vals[2] == "-"
                    });
            }
        };
    
    if (graph_pos) {
        // if we're given a graph_pos, we'll convert it into a path pos
        add_graph_pos(args::get(graph_pos));
    } else if (graph_pos_file) {
        std::ifstream gpos(args::get(graph_pos_file).c_str());
        std::string buffer;
        while (std::getline(gpos, buffer)) {
            add_graph_pos(buffer);
        }
    } else if (path_pos) {
        // if given a path pos, we convert it into a path pos in our reference set
        add_path_pos(args::get(path_pos));
    } else if (path_pos_file) {
        // if we're given a file of path positions, we'll convert them all
        std::ifstream refs(args::get(ref_path_file).c_str());
        std::string buffer;
        while (std::getline(refs, buffer)) {
            add_path_pos(buffer);
        }
    }
    // todo: bed files

    uint64_t search_radius = _search_radius ? args::get(_search_radius) : 10000;

    // make an hash set of our ref path ids for quicker lookup
    hash_set<uint64_t> ref_path_set;
    for (auto& path : ref_paths) {
        ref_path_set.insert(as_integer(path));
    }

    auto get_graph_pos =
        [&](const path_pos_t& pos) {
            auto path_end = graph.path_end(pos.path);
            uint64_t walked = 0;
            for (step_handle_t s = graph.path_begin(pos.path);
                 s != path_end; s = graph.get_next_step(s)) {
                handle_t h = graph.get_handle_of_step(s);
                walked += graph.get_length(h);
                if (walked > pos.offset) {
                    return make_pos_t(graph.get_id(h), graph.get_is_reverse(h), walked - pos.offset);
                }
            }
            assert(false);
            return make_pos_t(0, false, 0);
        };

    auto get_offset_in_path =
        [&](const path_handle_t& path, const step_handle_t& target) {
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
    
    // for each position that we want to look up
#pragma omp parallel for
    for (auto& pos : graph_positions) {
        // go to the graph
        // do a little BFS, bounded by our limit
        handle_t start_handle = graph.get_handle(id(pos), is_rev(pos));
        uint64_t max_dist = 0;
        bool found_hit = false;
        step_handle_t ref_hit;
        odgi::algorithms::bfs(graph,
            [&](const handle_t& h, const uint64_t& r, const uint64_t& l, const uint64_t& d) {
                max_dist = std::max(max_dist, l);
                bool got_hit = false;
                step_handle_t hit;
                graph.for_each_step_on_handle(
                    h, [&](const step_handle_t& s) {
                           auto p = graph.get_path_handle_of_step(s);
                           if (!got_hit && ref_path_set.count(as_integer(p))) {
                               got_hit = true;
                               hit = s;
                           }
                       });
                if (got_hit) {
                    ref_hit = hit;
                    found_hit = true;
                }
            },
            [](const handle_t& h) { return false; },
            [](const handle_t& l, const handle_t& h) { return false; },
            [&max_dist,&search_radius,&found_hit](void) {
                return max_dist > search_radius || found_hit;
            },
            { start_handle },
            { },
            true); // use bidirectional search
        // now, if we found our hit, print
        if (found_hit) {
            path_handle_t p = graph.get_path_handle_of_step(ref_hit);
            uint64_t path_offset = get_offset_in_path(p, ref_hit);
            bool ref_is_rev = false;
            std::cout << id(pos) << "," << offset(pos) << "," << (is_rev(pos) ? "-" : "+") << "\t"
                      << graph.get_path_name(p) << "," << path_offset << "," << (ref_is_rev ? "-" : "+") << std::endl;
        }
    }
#pragma omp parallel for
    for (auto& path_pos : path_positions) {
        auto pos = get_graph_pos(path_pos);
        handle_t start_handle = graph.get_handle(id(pos), is_rev(pos));
        uint64_t max_dist = 0;
        bool found_hit = false;
        step_handle_t ref_hit;
        odgi::algorithms::bfs(graph,
            [&](const handle_t& h, const uint64_t& r, const uint64_t& l, const uint64_t& d) {
                max_dist = std::max(max_dist, l);
                bool got_hit = false;
                step_handle_t hit;
                graph.for_each_step_on_handle(
                    h, [&](const step_handle_t& s) {
                           auto p = graph.get_path_handle_of_step(s);
                           if (!got_hit && ref_path_set.count(as_integer(p))) {
                               got_hit = true;
                               hit = s;
                           }
                       });
                if (got_hit) {
                    ref_hit = hit;
                    found_hit = true;
                }
            },
            [](const handle_t& h) { return false; },
            [](const handle_t& l, const handle_t& h) { return false; },
            [&max_dist,&search_radius,&found_hit](void) {
                return max_dist > search_radius || found_hit;
            },
            { start_handle },
            { },
            true); // use bidirectional search
        // now, if we found our hit, print
        if (found_hit) {
            path_handle_t p = graph.get_path_handle_of_step(ref_hit);
            uint64_t path_offset = get_offset_in_path(p, ref_hit);
            bool ref_is_rev = false;
            std::cout << graph.get_path_name(path_pos.path) << "," << path_pos.offset << "," << (path_pos.is_rev ? "-" : "+") << "\t"
                      << graph.get_path_name(p) << "," << path_offset << "," << (ref_is_rev ? "-" : "+") << std::endl;
        }
    }

    // go to the graph
    // run a bfs forward and backwards until we hit a reference path
    // if we hit multiple, choose the one that's first in our path list
    // todo: option to collect all mappings
    // output: write the input position, the reference position, and the distance and orientation to it

    
    return 0;
}

static Subcommand odgi_position("position", "coordinate projections between nodes and paths",
                                PIPELINE, 3, main_position);


}
