#include "subcommand.hpp"
#include "odgi.hpp"
#include "IntervalTree.h"
//#include "gfakluge.hpp"
#include "args.hxx"
#include "split.hpp"
//#include "io_helper.hpp"
#include "threads.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_stats(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi stats";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("metrics describing variation graphs and their path relationships");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the variation graph from this file", {'i', "idx"});
    args::Flag summarize(parser, "summarize", "summarize the graph properties and dimensions", {'S', "summarize"});
    args::Flag base_content(parser, "base-content", "describe the base content of the graph", {'b', "base-content"});
    args::Flag path_coverage(parser, "coverage", "provide a histogram of path coverage over bases in the graph", {'C', "coverage"});
    args::Flag path_setcov(parser, "setcov", "provide a histogram of coverage over unique sets of paths", {'V', "set-coverage"});
    args::Flag path_setcov_count(parser, "setcountcov", "provide a histogram of coverage over counts of unique paths", {'Q', "set-count-coverage"});
    args::Flag path_multicov(parser, "multicov", "provide a histogram of coverage over unique multisets of paths", {'M', "multi-coverage"});
    args::Flag path_multicov_count(parser, "multicountcov", "provide a histogram of coverage over counts of paths", {'L', "multi-count-coverage"});
    args::ValueFlag<std::string> path_bedmulticov(parser, "BED", "for each BED entry, provide a table of path coverage over unique multisets of paths in the graph. Each unique multiset of paths overlapping a given BED interval is described in terms of its length relative to the total interval, the number of path traversals, and unique paths involved in these traversals.", {'B', "bed-multicov"});
    args::ValueFlag<std::string> path_delim(parser, "CHAR", "the part of each path name before this delimiter is a group identifier, which when specified will cause stats to be collected in a group-wise rather than path-wise fashion", {'D', "delim"});
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
            graph.deserialize(std::cin);
        } else {
            ifstream f(infile.c_str());
            graph.deserialize(f);
            f.close();
        }
    }
    if (args::get(summarize)) {
        uint64_t length_in_bp = 0, node_count = 0, edge_count = 0, path_count = 0;
        graph.for_each_handle([&](const handle_t& h) {
                length_in_bp += graph.get_length(h);
                ++node_count;
            });
        graph.for_each_edge([&](const edge_t& e) {
                ++edge_count;
                return true;
            });
        graph.for_each_path_handle([&](const path_handle_t& p) {
                ++path_count;
            });
        std::cerr << "length:\t" << length_in_bp << std::endl;
        std::cerr << "nodes:\t" << node_count << std::endl;
        std::cerr << "edges:\t" << edge_count << std::endl;
        std::cerr << "paths:\t" << path_count << std::endl;
    }
    if (args::get(base_content)) {
        std::vector<uint64_t> chars(256);
        graph.for_each_handle([&](const handle_t& h) {
                std::string seq = graph.get_sequence(h);
                for (auto c : seq) {
                    ++chars[c];
                }
            });
        for (uint64_t i = 0; i < 256; ++i) {
            if (chars[i]) {
                std::cout << (char)i << "\t" << chars[i] << std::endl;
            }
        }
    }
    if (args::get(path_coverage)) {
        std::map<uint64_t, uint64_t> full_histogram;
        std::map<uint64_t, uint64_t> unique_histogram;
        graph.for_each_handle([&](const handle_t& h) {
                std::vector<uint64_t> paths_here;
                graph.for_each_step_on_handle(h, [&](const step_handle_t& occ) {
                        paths_here.push_back(as_integer(graph.get_path(occ)));
                    });
                std::sort(paths_here.begin(), paths_here.end());
                std::vector<uint64_t> unique_paths = paths_here;
                unique_paths.erase(std::unique(unique_paths.begin(), unique_paths.end()), unique_paths.end());
                full_histogram[paths_here.size()] += graph.get_length(h);
                unique_histogram[unique_paths.size()] += graph.get_length(h);
            });
        std::cout << "type\tcov\tN" << std::endl;
        for (auto& p : full_histogram) {
            std::cout << "full\t" << p.first << "\t" << p.second << std::endl;
        }
        for (auto& p : unique_histogram) {
            std::cout << "uniq\t" << p.first << "\t" << p.second << std::endl;
        }
    }

    bool using_delim = !args::get(path_delim).empty();
    char delim = '\0';
    ska::flat_hash_map<std::string, uint64_t> path_group_ids;
    std::vector<std::string> path_groups;
    if (using_delim) {
        delim = args::get(path_delim).at(0);
        uint64_t i = 0;
        graph.for_each_path_handle(
            [&](const path_handle_t& p) {
                std::string group_name = split(graph.get_path_name(p), delim)[0];
                auto f = path_group_ids.find(group_name);
                if (f == path_group_ids.end()) {
                    path_group_ids[group_name] = i++;
                    path_groups.push_back(group_name);
                }
            });
    }

    auto get_path_name
        = (using_delim ?
           (std::function<std::string(const uint64_t&)>)
           [&](const uint64_t& id) { return path_groups[id]; }
           :
           (std::function<std::string(const uint64_t&)>)
           [&](const uint64_t& id) { return graph.get_path_name(as_path_handle(id)); });

    auto get_path_id
        = (using_delim ?
           (std::function<uint64_t(const path_handle_t&)>)
           [&](const path_handle_t& p) {
               std::string group_name = split(graph.get_path_name(p), delim)[0];
               return path_group_ids[group_name];
           }
           :
           (std::function<uint64_t(const path_handle_t&)>)
           [&](const path_handle_t& p) {
               return as_integer(p);
           });

    if (args::get(path_setcov)) {
        uint64_t total_length = 0;
        std::map<std::set<uint64_t>, uint64_t> setcov;
        graph.for_each_handle(
            [&](const handle_t& h) {
                std::set<uint64_t> paths_here;
                graph.for_each_step_on_handle(
                    h,
                    [&](const step_handle_t& occ) {
                        paths_here.insert(get_path_id(graph.get_path(occ)));
                    });
                size_t l = graph.get_length(h);
#pragma omp critical (setcov)
                {
                    setcov[paths_here] += l;
                    total_length += l;
                }
            }, true);
        std::cout << "length\tgraph.frac\tn.paths\tpath.set" << std::endl;
        for (auto& p : setcov) {
            std::cout << p.second << "\t"
                      << (double)p.second/(double)total_length << "\t"
                      << p.first.size() << "\t";
            for (auto& i : p.first) {
                std::cout << get_path_name(i) << ",";
            }
            std::cout << std::endl;
        }
    }

    if (args::get(path_setcov_count)) {
        uint64_t total_length = 0;
        std::map<uint64_t, uint64_t> setcov_count;
        graph.for_each_handle(
            [&](const handle_t& h) {
                std::set<uint64_t> paths_here;
                graph.for_each_step_on_handle(
                    h,
                    [&](const step_handle_t& occ) {
                        paths_here.insert(get_path_id(graph.get_path(occ)));
                    });
                size_t l = graph.get_length(h);
#pragma omp critical (setcov)
                {
                    setcov_count[paths_here.size()] += l;
                    total_length += l;
                }
            }, true);
        std::cout << "length\tgraph.frac\tunique.path.count" << std::endl;
        for (auto& p : setcov_count) {
            std::cout << p.second << "\t"
                      << (double)p.second/(double)total_length
                      << "\t" << p.first << std::endl;
        }
    }

    if (args::get(path_multicov)) {
        uint64_t total_length = 0;
        std::map<std::vector<uint64_t>, uint64_t> multisetcov;
        graph.for_each_handle(
            [&](const handle_t& h) {
                std::vector<uint64_t> paths_here;
                graph.for_each_step_on_handle(
                    h,
                    [&](const step_handle_t& occ) {
                        paths_here.push_back(get_path_id(graph.get_path(occ)));
                    });
                std::sort(paths_here.begin(), paths_here.end());
                size_t l = graph.get_length(h);
#pragma omp critical (multisetcov)
                {
                    multisetcov[paths_here] += l;
                    total_length += l;
                }
            }, true);
        std::cout << "length\tgraph.frac\nn.paths\tpath.multiset" << std::endl;
        for (auto& p : multisetcov) {
            std::cout << p.second << "\t"
                      << (double)p.second/(double)total_length << "\t"
                      << p.first.size() << "\t";
            bool first = true;
            for (auto& i : p.first) {
                std::cout << (first ? (first=false, "") :",") << get_path_name(i);
            }
            std::cout << std::endl;
        }
    }

    if (args::get(path_multicov_count)) {
        uint64_t total_length = 0;
        std::map<uint64_t, uint64_t> multisetcov_count;
        graph.for_each_handle(
            [&](const handle_t& h) {
                uint64_t paths_here = 0;
                graph.for_each_step_on_handle(
                    h,
                    [&](const step_handle_t& occ) {
                        ++paths_here;
                    });
                size_t l = graph.get_length(h);
#pragma omp critical (multisetcov)
                {
                    multisetcov_count[paths_here] += l;
                    total_length += l;
                }
            }, true);
        std::cout << "length\tgraph.frac\tpath.step.count" << std::endl;
        for (auto& p : multisetcov_count) {
            std::cout << p.second << "\t"
                      << (double)p.second/(double)total_length << "\t"
                      << p.first << std::endl;
        }
    }

    if (!args::get(path_bedmulticov).empty()) {
        std::string line;
        typedef std::map<std::vector<uint64_t>, uint64_t> setcov_t;
        typedef IntervalTree<uint64_t, std::pair<std::string, setcov_t*> > itree_t;
        map<std::string, itree_t::interval_vector> intervals;
        auto& x = args::get(path_bedmulticov);
        std::ifstream bed_in(x);
        while (std::getline(bed_in, line)) {
            // BED is base-numbered, 0-origin, half-open.  This parse turns that
            // into base-numbered, 0-origin, fully-closed for internal use.  All
            // coordinates used internally should be in the latter, and coordinates
            // from the user in the former should be converted immediately to the
            // internal format.
            std::vector<string> fields = split(line, '\t');
            intervals[fields[0]].push_back(
                itree_t::interval(
                    std::stoul(fields[1]),
                    std::stoul(fields[2]),
                    make_pair(fields[3], new setcov_t())));
        }

        std::vector<std::string> path_names;
        path_names.reserve(intervals.size());
        for(auto const& i : intervals) {
            path_names.push_back(i.first);
        }

        // the header
        std::cout << "path.name" << "\t"
                  << "bed.name" << "\t"
                  << "bed.start" << "\t"
                  << "bed.stop" << "\t"
                  << "bed.len" << "\t"
                  << "path.set.state.len" << "\t"
                  << "path.set.frac" << "\t"
                  << "path.traversals" << "\t"
                  << "uniq.paths.in.state" << "\t"
                  << "path.multiset" << std::endl;

#pragma omp parallel for
        for (uint64_t k = 0; k < path_names.size(); ++k) {
            auto& path_name = path_names.at(k);
            auto& path_ivals = intervals[path_name];
            path_handle_t path = graph.get_path_handle(path_name);
            // build the intervals for each path we'll query
            itree_t itree(std::move(path_ivals)); //, 16, 1);
            uint64_t pos = 0;
            graph.for_each_step_in_path(graph.get_path_handle(x), [&](const step_handle_t& occ) {
                    std::vector<uint64_t> paths_here;
                    handle_t h = graph.get_handle_of_step(occ);
                    graph.for_each_step_on_handle(
                        h,
                        [&](const step_handle_t& occ) {
                            paths_here.push_back(get_path_id(graph.get_path(occ)));
                        });
                    uint64_t len = graph.get_length(h);
                    // check each position in the node
                    auto hits = itree.findOverlapping(pos, pos+len);
                    if (hits.size()) {
                        for (auto& h : hits) {
                            auto& q = *h.value.second;
                            // adjust length for overlap length
                            uint64_t ovlp = len - (h.start > pos ? h.start - pos : 0) - (h.stop < pos+len ? pos+len - h.stop : 0);
                            q[paths_here] += ovlp;
                        }
                    }
                    pos += len;
                });
            itree.visit_all([&](const itree_t::interval& ival) {
                    auto& name = ival.value.first;
                    auto& setcov = *ival.value.second;
                    for (auto& p : setcov) {
                        std::set<uint64_t> u(p.first.begin(), p.first.end());
#pragma omp critical (cout)
                        {
                            std::cout << path_name << "\t" << name << "\t" << ival.start << "\t" << ival.stop << "\t"
                                      << ival.stop - ival.start << "\t"
                                      << p.second << "\t"
                                      << (float)p.second/(ival.stop-ival.start) << "\t"
                                      << p.first.size() << "\t"
                                      << u.size() << "\t";
                            bool first = true;
                            for (auto& i : p.first) {
                                std::cout << (first ? (first=false, "") :",") << get_path_name(i);
                            }
                            std::cout << std::endl;
                        }
                    }
            });
            for (auto& ival : path_ivals) {
                delete ival.value.second;
            }
        }
    }

    return 0;
}

static Subcommand odgi_stats("stats", "describe the graph and its path relationships",
                              PIPELINE, 3, main_stats);


}
