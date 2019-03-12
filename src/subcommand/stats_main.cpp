#include "subcommand.hpp"
#include "graph.hpp"
//#include "gfakluge.hpp"
#include "args.hxx"
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
    args::Flag path_coverage(parser, "path_coverage", "provide a histogram of path coverage over bases in the graph", {'C', "path-coverage"});
    args::Flag path_setcov(parser, "path_setcov", "provide a histogram of coverage over unique sets of paths", {'V', "path-setcov"});
    args::Flag path_multicov(parser, "path_setcov", "provide a histogram of coverage over unique multisets of paths", {'M', "path-multicov"});
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

    return 0;
}

static Subcommand odgi_stats("stats", "extract statistics and properties of the graph",
                              PIPELINE, 3, main_stats);


}
