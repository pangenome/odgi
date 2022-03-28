#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/bubble.hpp"
#include "algorithms/stepindex.hpp"
#include "utils.hpp"
#include "split.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_bubble(int argc, char **argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc - 1; ++i) {
        argv[i] = argv[i + 1];
    }
    const std::string prog_name = "odgi bubble";
    argv[0] = (char *) prog_name.c_str();
    --argc
;
    args::ArgumentParser parser("Extract matrix of path pangenome coverage permutations for power law regression.");
    args::Group mandatory_opts(parser, "[ MANDATORY ARGUMENTS ]");
    args::ValueFlag<std::string> og_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::Group bubble_opts(parser, "[ Bubble Options ]");
    args::ValueFlag<std::string> _target_path(bubble_opts, "NAME", "Target (reference) path name.", {'p', "target-path"});
//    args::ValueFlag<std::string> _path_groups(bubble_opts, "FILE", "Group paths as described in two-column FILE, with columns path.name and group.name.",
//                                              {'p', "path-groups"});
//    args::Flag group_by_sample(bubble_opts, "bool", "Following PanSN naming (sample#hap#ctg), group by sample (1st field).", {'S', "group-by-sample"});
//    args::Flag group_by_haplotype(bubble_opts, "bool", "Following PanSN naming (sample#hap#ctg), group by haplotype (2nd field).", {'H', "group-by-haplotype"});
//    args::ValueFlag<std::string> _bed_targets(bubble_opts, "FILE", "BED file over path space of the graph, describing a subset of the graph to consider.",
//                                              {'b', "bed-targets"});
//    args::ValueFlag<uint64_t> _n_permutations(bubble_opts, "N", "Number of permutations to run.",
//                                             {'n', "n-permutations"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> nthreads(threading_opts, "N", "Number of threads to use for parallel operations.",
                                       {'t', "threads"});
    args::Group processing_info_opts(parser, "[ Processing Information ]");
    //args::Flag debug(processing_info_opts, "debug", "Print information about the process to stderr.", {'d', "debug"});
    args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi bubble.", {'h', "help"});
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
            << "[odgi::bubble] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
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
                utils::handle_gfa_odgi_input(infile, "bubble", args::get(progress), num_threads, graph);
            }
        }
    }

    graph.set_number_of_threads(num_threads);

    auto target_path = graph.get_path_handle(args::get(_target_path));
    
    algorithms::step_index_t step_index(graph, {target_path}, num_threads, true, 8);
    //step_index.save(step_index_out_file);

    //int last_depth = 0;
    auto handle_output = [&](const algorithms::bubble_t& bubble) {
        // 
    };

    algorithms::for_each_bubble(graph, step_index, target_path, handle_output);

    return 0;
}

static Subcommand odgi_bubble("bubble", "variant detection using linearized bubbles.",
                             PIPELINE, 3, main_bubble);


}
