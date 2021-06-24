#include "algorithms/untangle.hpp"
#include "args.hxx"
#include "odgi.hpp"
#include "subcommand.hpp"
#include "utils.hpp"
#include <omp.h>

namespace odgi {

using namespace odgi::subcommand;

int main_untangle(int argc, char **argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc - 1; ++i) {
        argv[i] = argv[i + 1];
    }
    const std::string prog_name = "odgi untangle";
    argv[0] = (char *)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("Resolve spurious inverting links.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> og_in_file(
        mandatory_opts, "FILE",
        "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually "
        "ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format "
        "requires additional time!",
        {'i', "idx"});
    args::Group untangling_opts(parser, "[ Untangling Options ]");
    args::ValueFlag<std::string> query(untangling_opts, "NAME", "Use this query path.",
                       {'q', "query-path"});
    args::ValueFlag<std::string> reference(untangling_opts, "NAME", "Use this reference path.",
                       {'r', "reference-path"});
    args::ValueFlag<uint64_t> merge_dist(untangling_opts, "N", "Merge segments shorter than this length into previous segments.",
                                         {'m', "merge-dist"});
    args::Group debugging_opts(parser, "[ Debugging Options ]");
    args::Flag make_self_dotplot(debugging_opts, "DOTPLOT", "Render a table showing the positional dotplot of the query against itself.",
                                 {'s', "self-dotplot"});
    args::Group threading(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> nthreads(
        threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
    args::Group processing_info_opts(parser, "[ Processing Information ]");
    args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.",
                        {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi untangle.",
                        {'h', "help"});

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
        std::cerr << "[odgi::untangle] error: please specify an input file from where to load the "
                     "graph via -i=[FILE], --idx=[FILE]."
                  << std::endl;
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
                utils::handle_gfa_odgi_input(infile, "untangle", args::get(progress), num_threads,
                                             graph);
            }
        }
    }

    if (make_self_dotplot) {
        algorithms::self_dotplot(graph, graph.get_path_handle(args::get(query)));
    } else {
        assert(query && reference);
        algorithms::untangle(graph,
                             { graph.get_path_handle(args::get(query)) },
                             { graph.get_path_handle(args::get(reference)) },
                             args::get(merge_dist),
                             num_threads);
    }

    return 0;
}

static Subcommand odgi_untangle("untangle", "Project paths into reference-relative BEDPE, to decompose paralogy relationships", PIPELINE, 3, main_untangle);

} // namespace odgi
