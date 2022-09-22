#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "split.hpp"
#include <omp.h>
#include "utils.hpp"
#include "algorithms/diffpriv.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_priv(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi priv";
    argv[0] = (char*)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("Differentially private transformation of graph path space.");
    args::Group mandatory_opts(parser, "[ MANDATORY ARGUMENTS ]");
    args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::ValueFlag<std::string> dg_out_file(mandatory_opts, "FILE", "Write the sorted dynamic succinct variation graph to this file. A file"
                                             " ending with *.og* is recommended.", {'o', "out"});
    args::Group mechanism_opts(parser, "[ Differential Privacy Mechanism Options ]");
    args::ValueFlag<double> input_epsilon(mechanism_opts, "e", "Epsilon for exponential mechanism.", {'e', "epsilon"});
    args::ValueFlag<double> target_depth(mechanism_opts, "DEPTH", "Sample until we have approximately this path depth.", {'d', "target-depth"});
    args::ValueFlag<uint64_t> input_min_hap_freq(mechanism_opts, "N", "Minimum frequency (count) of haplotype observation to emit. [default: 2]", {'c', "min-hap-freq"});
    args::ValueFlag<uint64_t> input_bp_limit(mechanism_opts, "bp", "Maximum sampled haplotype length. [default: 10000]", {'b', "bp-limit"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> threads(threading_opts, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
	args::Group processing_info_opts(parser, "[ Processing Information ]");
	args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi paths.", {'h', "help"});

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

    if (!dg_in_file) {
        std::cerr << "[odgi::priv] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (!dg_out_file || args::get(dg_out_file).empty()) {
        std::cerr << "[odgi::priv] error: please specify an output file to store the graph via -o=[FILE], --out=[FILE]." << std::endl;
        return 1;
    }

	const uint64_t num_threads = args::get(threads) ? args::get(threads) : 1;
    omp_set_num_threads(num_threads);

	graph_t graph;
    assert(argc > 0);
    std::string infile = args::get(dg_in_file);
    if (infile.size()) {
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
			utils::handle_gfa_odgi_input(infile, "paths", args::get(progress), num_threads, graph);
        }
    }

    double depth = target_depth ? args::get(target_depth) : 1;
    double epsilon = input_epsilon ? args::get(input_epsilon) : 0.001;
    uint64_t bp_limit = input_bp_limit ? args::get(input_bp_limit) : 10000;
    uint64_t min_haplotype_freq = input_min_hap_freq ? args::get(input_min_hap_freq) : 2;
    if (min_haplotype_freq < 2) {
        std::cerr << "[odgi::priv] ERROR: setting -c, --min-hap-freq to less than 2 is not supported due to singularities that arise in the exponential mechanism." << std::endl;
        return 1;
    }

    graph_t priv;

    algorithms::diff_priv(graph, priv, epsilon, depth, min_haplotype_freq, bp_limit, num_threads);

    const std::string outfile = args::get(dg_out_file);
    if (outfile == "-") {
        graph.serialize(std::cout);
    } else {
        ofstream f(outfile.c_str());
        graph.serialize(f);
        f.close();
    }

    return 0;
}

static Subcommand odgi_priv("priv", "Make a differentially-private transformation of graph paths.",
                            PIPELINE, 3, main_priv);


}
