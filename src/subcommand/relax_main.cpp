#include "subcommand.hpp"
#include <iostream>
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/layout.hpp"
#include <numeric>
#include "progress.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_relax(int argc, char **argv) {

    // trick argument parser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc - 1; ++i) {
        argv[i] = argv[i + 1];
    }
    std::string prog_name = "odgi relax";
    argv[0] = (char *) prog_name.c_str();
    --argc;

    args::ArgumentParser parser(
        "relax a graph by converting global large structural variants into large local bubbles");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> bed_in_file(parser, "FILE", "read the sorted BED file from this file", {'B', "bed-in"});
    args::Flag progress(parser, "progress", "display progress", {'P', "progress"});
    args::ValueFlag<uint64_t> nthreads(parser, "N", "number of threads to use for parallel phases", {'t', "threads"});

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

    uint64_t thread_count = 1;
    if (args::get(nthreads)) {
        omp_set_num_threads(args::get(nthreads));
        thread_count = args::get(nthreads);
    }

    if (!dg_in_file) {
        std::cerr
                << "[odgi relax] error: Please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
                << std::endl;
        return 1;
    }

    if (!bed_in_file) {
        std::cerr
                << "[odgi relax] error: Please specify an input file from where to read the BED records via -B=[FILE], --bed-in=[FILE]."
                << std::endl;
        return 1;
    }

    graph_t graph;
    assert(argc > 0);
    std::string infile = args::get(dg_in_file);
    if (infile.size()) {
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
            ifstream f(infile.c_str());
            graph.deserialize(f);
            f.close();
        }
    }

    // TODO read BED line by line and put each of them into a min_bed_record_t struct:
    //        std::string chrom;
    //        uint64_t chromStart;
    //        uint64_t chromEnd;
    //        double path_layout_nuc_dist_ratio;
    // The struct will be made in bed_records.hpp
    // TODO calculate median path_layout_nuc_dist_ratio
    // TODO only keep min_bed_record_t structs which are double -m/--min-median-factor above the median
    // TODO only keep the -r/--relax-num OR -e/--relax-percentage min_bed_record_t structs
    // TODO add these as a vector<min_bed_record_t> to a flat hash map for each path
    // the flat hash map is the result of the function, which will be stored here

    // TODO iterate over the entries of the hash map, in parallel, if possible, and relax
    // 1. iterate over all steps, if we have a hit within the given range, note down the handles and start and end step, careful about orientation!
    // 2. remove the steps from the graph, if possible also the nodes and edges
    // 3. add new nodes and steps using std::pair<step_handle_t, step_handle_t> graph_t::rewrite_segment(const step_handle_t& segment_begin,
    //                                                                 const step_handle_t& segment_end,
    //                                                                 const std::vector<handle_t>& new_segment)

    // TODO clean up the graph, removing zero coverage nodes and edges

    // TODO write out the graph

    return 0;
}

static Subcommand odgi_relax("relax", "relax a graph by converting global large structural variants into large local bubbles",
                            PIPELINE, 3, main_relax);


}
