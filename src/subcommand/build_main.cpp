#include "subcommand.hpp"
#include "graph.hpp"
#include "gfakluge.hpp"
#include "args.hxx"
//#include "io_helper.hpp"

namespace dg {

using namespace dg::subcommand;

int main_build(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "dg build";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("construct a dynamic succinct variation graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> gfa_file(parser, "FILE", "construct the graph from this GFA input file", {'g', "gfa"});
    args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the index in this file", {'o', "out"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the index from this file", {'i', "idx"});
    //args::ValueFlag<std::string> seqs(parser, "FILE", "the sequences used to generate the alignments", {'s', "seqs"});
    //args::ValueFlag<std::string> base(parser, "FILE", "build graph using this basename", {'b', "base"});
    //args::ValueFlag<uint64_t> num_threads(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    //args::ValueFlag<uint64_t> repeat_max(parser, "N", "limit transitive closure to include no more than N copies of a given input base", {'r', "repeat-max"});
    //args::ValueFlag<uint64_t> aln_keep_n_longest(parser, "N", "keep up to the N-longest alignments overlapping each query position", {'k', "aln-keep-n-longest"});
    //args::ValueFlag<uint64_t> aln_min_length(parser, "N", "ignore alignments shorter than this", {'m', "aln-min-length"});
    args::Flag to_gfa(parser, "to_gfa", "write the graph to stdout in GFA format", {'G', "to-gfa"});
    args::Flag summarize(parser, "summarize", "summarize the graph properties and dimensions", {'S', "summarize"});
    args::Flag debug(parser, "debug", "enable debugging", {'d', "debug"});
    args::Flag progress(parser, "progress", "show progress updates", {'p', "progress"});
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

    /*
    size_t n_threads = args::get(num_threads);
    if (n_threads) {
        omp_set_num_threads(args::get(num_threads));
    } else {
        omp_set_num_threads(1);
    }
    */
    graph_t graph;
    
    //make_graph();
    assert(argc > 0);
    std::string gfa_filename = args::get(gfa_file);
    if (gfa_filename.size()) {
        char* filename = (char*)gfa_filename.c_str();
        //std::cerr << "filename is " << filename << std::endl;
        gfak::GFAKluge gg;
        //double version = gg.detect_version_from_file(filename);
        //std::cerr << version << " be version" << std::endl;
        //assert(version == 1.0);
        /*
          uint64_t num_nodes = 0;
          gg.for_each_sequence_line_in_file(filename, [&](gfak::sequence_elem s) {
          ++num_nodes;
          });
          graph_t graph(num_nodes+1); // include delimiter
        */
        uint64_t i = 0;
        gg.for_each_sequence_line_in_file(filename, [&](gfak::sequence_elem s) {
                uint64_t id = stol(s.name);
                graph.create_handle(s.sequence, id);
                if (args::get(progress)) {
                    if (i % 1000 == 0) std::cerr << "node " << i << "\r";
                    ++i;
                }
            });
        if (args::get(progress)) {
            i = 0; std::cerr << std::endl;
        }
        gg.for_each_edge_line_in_file(filename, [&](gfak::edge_elem e) {
                if (e.source_name.empty()) return;
                handle_t a = graph.get_handle(stol(e.source_name), !e.source_orientation_forward);
                handle_t b = graph.get_handle(stol(e.sink_name), !e.sink_orientation_forward);
                graph.create_edge(a, b);
                if (args::get(progress)) {
                    if (i % 1000 == 0) std::cerr << "edge " << i << "\r";
                    ++i;
                }
            });
        if (args::get(progress)) {
            i = 0; std::cerr << std::endl;
        }
        gg.for_each_path_element_in_file(filename, [&](const std::string& path_name, const std::string& node_id, bool is_rev, const std::string& cigar) {
                path_handle_t path;
                if (!graph.has_path(path_name)) {
                    path = graph.create_path_handle(path_name);
                    if (args::get(progress) && i != 0) {
                        std::cerr << std::endl;
                    }
                    i = 0;
                } else {
                    path = graph.get_path_handle(path_name);
                }
                handle_t occ = graph.get_handle(stol(node_id), is_rev);
                graph.append_occurrence(path, occ);
                if (args::get(progress)) {
                    if (i % 1000 == 0) std::cerr << "path " << path_name << " " << i << "\r";
                    ++i;
                }
                // ignores overlaps
            });
    }
    std::string infile = args::get(dg_in_file);
    if (infile.size()) {
        ifstream f(infile.c_str());
        graph.load(f);
        f.close();
    }
    if (args::get(progress)) {
        std::cerr << std::endl;
    }
    // here we should measure memory usage etc.
    if (args::get(debug)) {
        graph.display();
    }
    if (args::get(to_gfa)) {
        graph.to_gfa(std::cout);
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
    std::string outfile = args::get(dg_out_file);
    if (outfile.size()) {
        ofstream f(outfile.c_str());
        graph.serialize(f);
        f.close();
    }
    //if (args::get(
    return 0;
}

static Subcommand dg_build("build", "build dynamic succinct variation graph",
                           PIPELINE, 3, main_build);


}
