#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/normalize.hpp"
#include "utils.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_normalize(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    const std::string prog_name = "odgi normalize";
    argv[0] = (char*)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("Compact unitigs and simplify redundant furcations.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> og_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::ValueFlag<std::string> og_out_file(mandatory_opts, "FILE", "Write the normalized dynamic succinct variation graph in ODGI format to this file. A"
                                                                     " file ending with *.og* is recommended.", {'o', "out"});
    args::Group norm_opts(parser, "[ Normalize Options ]");
    args::ValueFlag<uint64_t> max_iterations(norm_opts, "N", "Iterate the normalization up to N many times (default: 10).", {'I', "max-iterations"});
	args::Group threading(parser, "[ Threading ]");
	args::ValueFlag<uint64_t> nthreads(threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
    args::Group process_info_opts(parser, "[ Processing Information ]");
    args::Flag debug(process_info_opts, "debug", "Print information about the normalization process to stdout.", {'d', "debug"});
	args::Flag progress(process_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
	args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi normalize.", {'h', "help"});

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
        std::cerr << "[odgi::normalize] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (!og_out_file) {
        std::cerr << "[odgi::normalize] error: please specify an output file to where to store the normalized graph via -o=[FILE], --out=[FILE]." << std::endl;
        return 1;
    }

	const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

	graph_t graph;
    assert(argc > 0);
    {
        const std::string infile = args::get(og_in_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
				utils::handle_gfa_odgi_input(infile, "normalize", args::get(progress), num_threads, graph);
            }
        }
    }

    algorithms::normalize(graph, args::get(max_iterations) ? args::get(max_iterations) : 10, args::get(debug));

    {
        const std::string outfile = args::get(og_out_file);
        if (!outfile.empty()) {
            if (outfile == "-") {
                graph.serialize(std::cout);
            } else {
                ofstream f(outfile.c_str());
                graph.serialize(f);
                f.close();
            }
        }
    }

    return 0;
}

static Subcommand odgi_normalize("normalize", "Compact unitigs and simplify redundant furcations.",
                                 PIPELINE, 3, main_normalize);


}
