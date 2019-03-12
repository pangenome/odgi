#include "subcommand.hpp"
#include "graph.hpp"
#include "IntervalTree.h"
//#include "gfakluge.hpp"
#include "args.hxx"
#include "split.hpp"
//#include "io_helper.hpp"

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
    
    args::ArgumentParser parser("metrics describing variation graphs");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    //args::ValueFlag<std::string> gfa_file(parser, "FILE", "construct the graph from this GFA input file", {'g', "gfa"});
    //args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the index in this file", {'o', "out"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the index from this file", {'i', "idx"});
    //args::ValueFlag<std::string> seqs(parser, "FILE", "the sequences used to generate the alignments", {'s', "seqs"});
    //args::ValueFlag<std::string> base(parser, "FILE", "build graph using this basename", {'b', "base"});
    //args::ValueFlag<uint64_t> num_threads(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    //args::ValueFlag<uint64_t> repeat_max(parser, "N", "limit transitive closure to include no more than N copies of a given input base", {'r', "repeat-max"});
    //args::ValueFlag<uint64_t> aln_keep_n_longest(parser, "N", "keep up to the N-longest alignments overlapping each query position", {'k', "aln-keep-n-longest"});
    //args::ValueFlag<uint64_t> aln_min_length(parser, "N", "ignore alignments shorter than this", {'m', "aln-min-length"});
    //args::Flag to_gfa(parser, "to_gfa", "write the graph to stdout in GFA format", {'G', "to-gfa"});
    args::Flag summarize(parser, "summarize", "summarize the graph properties and dimensions", {'S', "summarize"});
    args::Flag path_coverage(parser, "coverage", "provide a histogram of path coverage over bases in the graph", {'C', "coverage"});
    args::Flag path_setcov(parser, "setcov", "provide a histogram of coverage over unique sets of paths", {'V', "set-coverage"});
    args::Flag path_multicov(parser, "multicov", "provide a histogram of coverage over unique multisets of paths", {'M', "multi-coverage"});
    args::ValueFlag<std::string> path_bedmulticov(parser, "BED", "for each BED entry, provide a histogram of coverage over unique multisets of paths", {'B', "bed-multicov"});
    //args::Flag debug(parser, "debug", "enable debugging", {'d', "debug"});
    //args::Flag progress(parser, "progress", "show progress updates", {'p', "progress"});
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

    graph_t graph;
    assert(argc > 0);
    std::string infile = args::get(dg_in_file);
    if (infile.size()) {
        ifstream f(infile.c_str());
        graph.load(f);
        f.close();
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
    if (args::get(path_coverage)) {
        std::map<uint64_t, uint64_t> full_histogram;
        std::map<uint64_t, uint64_t> unique_histogram;
        graph.for_each_handle([&](const handle_t& h) {
                std::vector<uint64_t> paths_here;
                graph.for_each_occurrence_on_handle(h, [&](const occurrence_handle_t& occ) {
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
    if (args::get(path_setcov)) {
        std::map<std::set<uint64_t>, uint64_t> setcov;
        graph.for_each_handle([&](const handle_t& h) {
                std::set<uint64_t> paths_here;
                graph.for_each_occurrence_on_handle(h, [&](const occurrence_handle_t& occ) {
                        paths_here.insert(as_integer(graph.get_path(occ)));
                    });
                setcov[paths_here] += graph.get_length(h);
            });
        std::cout << "cov\tsets" << std::endl;
        for (auto& p : setcov) {
            std::cout << p.second << "\t";
            for (auto& i : p.first) {
                std::cout << graph.get_path_name(as_path_handle(i)) << ",";
            }
            std::cerr << std::endl;
        }
    }
    if (args::get(path_multicov)) {
        std::map<std::vector<uint64_t>, uint64_t> setcov;
        graph.for_each_handle([&](const handle_t& h) {
                std::vector<uint64_t> paths_here;
                graph.for_each_occurrence_on_handle(h, [&](const occurrence_handle_t& occ) {
                        paths_here.push_back(as_integer(graph.get_path(occ)));
                    });
                std::sort(paths_here.begin(), paths_here.end());
                setcov[paths_here] += graph.get_length(h);
            });
        std::cout << "cov\tsets" << std::endl;
        for (auto& p : setcov) {
            std::cout << p.second << "\t";
            for (auto& i : p.first) {
                std::cout << graph.get_path_name(as_path_handle(i)) << ",";
            }
            std::cout << std::endl;
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

        for (auto& p : intervals) {
            auto& path_name = p.first;
            path_handle_t path = graph.get_path_handle(path_name);
            // build the intervals for each path we'll query
            itree_t itree(std::move(p.second)); //, 16, 1);
            uint64_t pos = 0;
            graph.for_each_occurrence_in_path(graph.get_path_handle(x), [&](const occurrence_handle_t& occ) {
                    std::vector<uint64_t> paths_here;
                    handle_t h = graph.get_occurrence(occ);
                    graph.for_each_occurrence_on_handle(h, [&](const occurrence_handle_t& occ) {
                            paths_here.push_back(as_integer(graph.get_path(occ)));
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
                        std::cout << path_name << "\t" << name << "\t" << ival.start << "\t" << ival.stop << "\t"
                                  << p.second << "\t"
                                  << u.size() << "\t";
                        for (auto& i : p.first) {
                            std::cout << graph.get_path_name(as_path_handle(i)) << ",";
                        }
                        std::cout << std::endl;
                    }
            });
            for (auto& ival : p.second) {
                delete ival.value.second;
            }
        }
    }

    return 0;
}

static Subcommand odgi_stats("stats", "extract statistics and properties of the graph",
                              PIPELINE, 3, main_stats);


}
