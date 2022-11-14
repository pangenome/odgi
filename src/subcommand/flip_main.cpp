#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/flip.hpp"
#include "utils.hpp"
#include "split.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_flip(int argc, char **argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc - 1; ++i) {
        argv[i] = argv[i + 1];
    }
    const std::string prog_name = "odgi flip";
    argv[0] = (char *) prog_name.c_str();
    --argc
;
    args::ArgumentParser parser("Flip (reverse complement) paths to match the graph.");
    args::Group mandatory_opts(parser, "[ MANDATORY ARGUMENTS ]");
    args::ValueFlag<std::string> og_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. GFA is also supported.", {'i', "idx"});
    args::ValueFlag<std::string> og_out_file(mandatory_opts, "FILE", "Write the sorted dynamic succinct variation graph to this file (e.g. *.og*).", {'o', "out"});
    args::Group flip_opts(parser, "[ Flip Options ]");
    args::ValueFlag<std::string> _no_flips(flip_opts, "FILE", "Don't flip paths listed one per line in FILE.", {'n', "no-flips"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> nthreads(threading_opts, "N", "Number of threads to use for parallel operations.",
                                       {'t', "threads"});
    args::Group processing_info_opts(parser, "[ Processing Information ]");
    //args::Flag debug(processing_info_opts, "debug", "Print information about the process to stderr.", {'d', "debug"});
    args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi flip.", {'h', "help"});
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
    if (argc == 1) {
        std::cout << parser;
        return 1;
    }

    if (!og_in_file) {
        std::cerr
            << "[odgi::flip] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
            << std::endl;
        return 1;
    }

    if (!og_out_file || args::get(og_out_file).empty()) {
        std::cerr << "[odgi::flip] error: please specify an output file to store the graph via -o=[FILE], --out=[FILE]." << std::endl;
        return 1;
    }

    const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;
    omp_set_num_threads(num_threads);

    graph_t graph;
    assert(argc > 0);
    {
        const std::string infile = args::get(og_in_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
                utils::handle_gfa_odgi_input(infile, "flip", args::get(progress), num_threads, graph);
            }
        }
    }

    graph.set_number_of_threads(num_threads);

    std::vector<path_handle_t> no_flips;
    if (_no_flips) {
        no_flips = utils::load_paths(args::get(_no_flips), graph, "flip");
    }

    graph_t into;
    algorithms::flip_paths(graph, into, no_flips);

    const std::string outfile = args::get(og_out_file);
    if (outfile == "-") {
        into.serialize(std::cout);
    } else {
        ofstream f(outfile.c_str());
        into.serialize(f);
        f.close();
    }

    return 0;
}

static Subcommand odgi_flip("flip", "Flip path orientations to match the graph.",
                              PIPELINE, 3, main_flip);


}
