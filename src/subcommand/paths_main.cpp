#include "subcommand.hpp"
#include "graph.hpp"
//#include "IntervalTree.h"
//#include "gfakluge.hpp"
#include "args.hxx"
#include "split.hpp"
//#include "io_helper.hpp"
#include "position.hpp"
#include "threads.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_paths(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi paths";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("embedded path interrogation");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::Flag list_names(parser, "list-names", "list the paths in the graph", {'L', "list-paths"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the index from this file", {'i', "idx"});
    args::ValueFlag<std::string> overlaps_file(parser, "FILE", "Each line in (tab-delimited) FILE lists a grouping and a path. For each group we will provide pairwise overlap statistics for each pairing.", {'O', "overlaps"});
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

    if (args::get(threads)) {
        omp_set_num_threads(args::get(threads));
    } else {
        omp_set_num_threads(1);
    }

    graph_t graph;
    assert(argc > 0);
    std::string infile = args::get(dg_in_file);
    if (infile.size()) {
        if (infile == "-") {
            graph.load(std::cin);
        } else {
            ifstream f(infile.c_str());
            graph.load(f);
            f.close();
        }
    }

    //args::Flag list_names(parser, "list-names", "list the paths in the graph", {'L', "list-paths"});
    if (args::get(list_names)) {
        graph.for_each_path_handle([&](const path_handle_t& p) {
                std::cout << graph.get_path_name(p) << std::endl;
            });
    }
    
    if (!args::get(overlaps_file).empty()) {
        std::string line;
        ska::flat_hash_map<std::string, std::vector<std::string> > path_sets;
        ska::flat_hash_map<std::string, std::vector<pos_t> > path_decomposition;
        auto& x = args::get(overlaps_file);
        std::ifstream overlaps_in(x);
        while (std::getline(overlaps_in, line)) {
            // This file should contain tab-delimited lists of path names, one per line
            // the first field is an identifier
            std::vector<string> fields = split(line, '\t');
            assert(fields.size()==2);
            auto& name = fields[0];
            auto& pathset = path_sets[name];
            pathset.push_back(fields[1]);
            path_decomposition[fields[1]] = {};
        }

        std::vector<std::string> path_names;
        for (auto& p : path_decomposition) {
            path_decomposition[p.first] = {};
            path_names.push_back(p.first);
        }

        // the header
        std::cout << "group.name" << "\t"
                  << "query" << "\t"
                  << "target" << "\t"
                  << "overlap" << "\t"
                  << "overlap.frac" << std::endl;

#pragma omp parallel for
        for (uint64_t k = 0; k < path_names.size(); ++k) {
            auto& path_name = path_names.at(k);
            auto& decomposition = path_decomposition[path_name];
            // walk the path, adding each position to the decomposition
            path_handle_t path = graph.get_path_handle(path_name);
            uint64_t pos = 0;
            graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
                    handle_t h = graph.get_handle_of_step(occ);
                    nid_t id = graph.get_id(h);
                    uint64_t len = graph.get_length(h);
                    for (uint64_t i = 0; i < len; ++i) {
                        decomposition.push_back(make_pos_t(id, i, graph.get_is_reverse(h)));
                    }
                });
        }

        for (auto& s : path_sets) {
            auto& group_name = s.first;
            auto& paths = s.second;
            for (uint64_t i = 0; i < paths.size(); ++i) {
                for (uint64_t j = i+1; j < paths.size(); ++j) {
                    auto& v1 = path_decomposition[paths[i]];
                    auto& v2 = path_decomposition[paths[j]];
                    std::vector<pos_t> v3;
                    std::sort(v1.begin(), v1.end());
                    std::sort(v2.begin(), v2.end());
                    std::set_intersection(v1.begin(),v1.end(),
                                          v2.begin(),v2.end(),
                                          back_inserter(v3));
                    //ska::flat_hash_map<std::string, std::vector<pos_t> > path_decomposition;                    
                    std::cout << group_name << "\t" << paths[i] << "\t" << paths[j]
                              << "\t" << v3.size()
                              << "\t" << (float)v3.size()/((float)(v1.size()+v2.size())/2) << std::endl;
                }
            }
        }
    }

    return 0;
}

static Subcommand odgi_paths("paths", "interrogation and manipulation of paths",
                              PIPELINE, 3, main_paths);


}
