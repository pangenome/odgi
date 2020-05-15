#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "algorithms/linear_index.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_flatten(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi flatten";
    argv[0] = (char*)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("generate linearizations of the graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> odgi_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> fasta_out_file(parser, "FILE", "write concatenated node sequences in FASTA format to FILE", {'f', "fasta"});
    args::ValueFlag<std::string> fasta_seq_name(parser, "FILE", "name to use for the concatenated graph sequence (default: input file name)", {'n', "name-seq"});
    args::ValueFlag<std::string> bed_out_file(parser, "FILE", "write a BED format FILE describing the mapping between graph paths and the linearized FASTA sequence", {'b', "bed"});

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

    if (!odgi_in_file) {
        std::cerr << "[odgi flatten] error: Please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    graph_t graph;
    assert(argc > 0);
    std::string infile = args::get(odgi_in_file);
    if (infile.size()) {
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
            ifstream f(infile.c_str());
            graph.deserialize(f);
            f.close();
        }
    }

    if (args::get(fasta_out_file).empty() && args::get(bed_out_file).empty()) {
        std::cerr << "[odgi flatten] Error: a FASTA or BED output must be specified" << std::endl;
        return 2;
    }

    // graph linearization with handle to position mapping
    algorithms::linear_index_t linear(graph);

    std::string fasta_name = !args::get(fasta_seq_name).empty() ? args::get(fasta_seq_name) : args::get(odgi_in_file);
    std::string fasta_out = args::get(fasta_out_file);
    uint64_t fasta_line_width = 80;
    auto write_fasta
        = [&](ostream& out) {
              out << ">" << fasta_name << std::endl;
              for (uint64_t i = 0; i < linear.graph_seq.size(); i+=fasta_line_width) {
                  out << linear.graph_seq.substr(i, fasta_line_width) << std::endl;
              }
          };
    if (fasta_out.size()) {
        if (fasta_out == "-") {
            write_fasta(std::cout);
        } else {
            ofstream f(fasta_out.c_str());
            write_fasta(f);
        }
    }

    if (!args::get(bed_out_file).empty()) {
        auto write_bed_line
            = [&](ostream& out, const std::string& path_name, const uint64_t& start,
                  const uint64_t& end, bool is_rev, const uint64_t& rank) {
                  out << fasta_name << "\t" << start << "\t" << end << "\t"
                      << path_name << "\t" << (is_rev ? "-" : "+") << "\t" << rank << std::endl;
              };
        std::string bed_out = args::get(bed_out_file);
        ofstream b;
        bool bed_stdout = true;
        if (bed_out != "-") {
            b.open(bed_out.c_str());
            bed_stdout = false;
        }
        graph.for_each_path_handle(
            [&](const path_handle_t& p) {
                std::string path_name = graph.get_path_name(p);
                uint64_t rank = 0;
                graph.for_each_step_in_path(
                    p,
                    [&](const step_handle_t& s) {
                        handle_t h = graph.get_handle_of_step(s);
                        uint64_t start = linear.position_of_handle(h);
                        uint64_t end = start + graph.get_length(h);
                        bool is_rev = graph.get_is_reverse(h);
                        if (bed_stdout) {
                            write_bed_line(std::cout, path_name, start, end, is_rev, rank++);
                        } else {
                            write_bed_line(b, path_name, start, end, is_rev, rank++);
                        }
                    });
            });
    }
    
    return 0;
}

static Subcommand odgi_flatten("flatten", "project the graph sequence and paths into FASTA and BED",
                               PIPELINE, 3, main_flatten);


}
