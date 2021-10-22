#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/linear_index.hpp"
#include "utils.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_flatten(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    const std::string prog_name = "odgi flatten";
    argv[0] = (char*)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("Generate linearizations of a graph.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> odgi_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::Group output_opts(parser, "[ Output Options ]");
    args::ValueFlag<std::string> fasta_out_file(output_opts, "FILE", "Write the concatenated node sequences (also known as pangenome sequence) in FASTA format to FILE.", {'f', "fasta"});
    args::ValueFlag<std::string> fasta_seq_name(output_opts, "FILE", "The name to use for the concatenated graph sequence (default: input file name which was specified via -i, --idx=[FILE]).", {'n', "name-seq"});
    args::ValueFlag<std::string> bed_out_file(output_opts, "FILE", "Write the mapping between graph paths and the linearized FASTA sequence in BED format to FILE.", {'b', "bed"});
	args::Group threading(parser, "[ Threading ]");
	args::ValueFlag<uint64_t> nthreads(threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
	args::Group processing_info_opts(parser, "[ Processing Information ]");
	args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi flatten.", {'h', "help"});

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
        std::cerr << "[odgi::flatten] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (args::get(fasta_out_file).empty() && args::get(bed_out_file).empty()) {
        std::cerr << "[odgi::flatten] error: a FASTA or BED output must be specified" << std::endl;
        return 2;
    }

	const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

	graph_t graph;
    assert(argc > 0);
    {
        const std::string infile = args::get(odgi_in_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
				utils::handle_gfa_odgi_input(infile, "flatten", args::get(progress), num_threads, graph);;
            }
        }
    }

    // graph linearization with handle to position mapping
    algorithms::linear_index_t linear(graph);

    const std::string fasta_name = !args::get(fasta_seq_name).empty() ? args::get(fasta_seq_name) : args::get(odgi_in_file);

    {
        const std::string fasta_out = args::get(fasta_out_file);
        if (!fasta_out.empty()) {
            const uint8_t fasta_line_width = 80;

            auto write_fasta = [&](ostream& out) {
                out << ">" << fasta_name << std::endl;
                for (uint64_t i = 0; i < linear.graph_seq.size(); i+=fasta_line_width) {
                    out << linear.graph_seq.substr(i, fasta_line_width) << std::endl;
                }
            };

            if (fasta_out == "-") {
                write_fasta(std::cout);
            } else {
                ofstream f(fasta_out.c_str());
                write_fasta(f);
            }
        }
    }

    if (!args::get(bed_out_file).empty()) {
        auto write_bed_line = [&](ostream& out, const std::string& path_name, const uint64_t& start,
                const uint64_t& end, bool is_rev, const uint64_t& rank) {
            out << fasta_name << "\t" << start << "\t" << end << "\t"
            << path_name << "\t" << (is_rev ? "-" : "+") << "\t" << rank << std::endl;
        };

        const std::string bed_out = args::get(bed_out_file);
        ofstream b;
        bool bed_stdout = true;
        if (bed_out != "-") {
            bed_stdout = false;
            b.open(bed_out.c_str());
        }

        (bed_stdout ? std::cout : b) << "#name\tstart\tend\tpath.name\tstrand\tstep.rank\n";

        graph.for_each_path_handle([&](const path_handle_t& p) {
            const std::string path_name = graph.get_path_name(p);
            uint64_t rank = 0;
            graph.for_each_step_in_path(p, [&](const step_handle_t& s) {
                const handle_t h = graph.get_handle_of_step(s);
                const uint64_t start = linear.position_of_handle(h);
                const uint64_t end = start + graph.get_length(h);
                const bool is_rev = graph.get_is_reverse(h);
                write_bed_line(bed_stdout ? std::cout : b, path_name, start, end, is_rev, rank++);
            });
        });

        if (!bed_stdout) {
            b.close();
        }
    }
    
    return 0;
}

static Subcommand odgi_flatten("flatten", "Generate linearizations of a graph.",
                               PIPELINE, 3, main_flatten);


}
