#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/procbed.hpp"
#include "utils.hpp"
#include "split.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_procbed(int argc, char **argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc - 1; ++i) {
        argv[i] = argv[i + 1];
    }
    const std::string prog_name = "odgi procbed";
    argv[0] = (char *) prog_name.c_str();
    --argc
;
    args::ArgumentParser parser("Pprocrustes-BED: Intersect and adjust BED interval into PanSN-defined path subranges. Lift BED files into graphs produced by odgi extract. Uses path range information in the path names.");
    args::Group mandatory_opts(parser, "[ MANDATORY ARGUMENTS ]");
    args::ValueFlag<std::string> og_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::Group procbed_opts(parser, "[ Procbed Options ]");
    args::ValueFlag<std::string> _bed_targets(procbed_opts, "FILE", "BED file over path space of the full graph from which this subgraph was obtained. Using path range information in the path names, overlapping records will be rewritten to fit into the coordinate space of the subgraph.",
                                              {'b', "bed-targets"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> nthreads(threading_opts, "N", "Number of threads to use for parallel operations.",
                                       {'t', "threads"});
    args::Group processing_info_opts(parser, "[ Processing Information ]");
    //args::Flag debug(processing_info_opts, "debug", "Print information about the process to stderr.", {'d', "debug"});
    args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi procbed.", {'h', "help"});
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
            << "[odgi::procbed] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
            << std::endl;
        return 1;
    }

    std::string bed_targets;
    if (_bed_targets) {
        bed_targets = args::get(_bed_targets);
    } else {
        std::cerr
            << "[odgi::procbed] error: please specify an input BED file to adjust -b=[FILE], --bed-targets=[FILE]."
            << std::endl;
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
                utils::handle_gfa_odgi_input(infile, "procbed", args::get(progress), num_threads, graph);
            }
        }
    }

    graph.set_number_of_threads(num_threads);

    algorithms::adjust_ranges(graph, bed_targets);

    return 0;
}

static Subcommand odgi_procbed("procbed", "Procrustes-BED: adjust BED to match subpaths in graph.",
                              PIPELINE, 3, main_procbed);


}
