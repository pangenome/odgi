#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/inject.hpp"
#include "utils.hpp"
#include "split.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_inject(int argc, char **argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc - 1; ++i) {
        argv[i] = argv[i + 1];
    }
    const std::string prog_name = "odgi inject";
    argv[0] = (char *) prog_name.c_str();
    --argc
;
    args::ArgumentParser parser("Inject BED interval ranges as paths in the graph.");
    args::Group mandatory_opts(parser, "[ MANDATORY ARGUMENTS ]");
    args::ValueFlag<std::string> og_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::ValueFlag<std::string> og_out_file(mandatory_opts, "FILE", "Write the sorted dynamic succinct variation graph to this file. A file ending with *.og* is recommended.", {'o', "out"});
    args::Group inject_opts(parser, "[ Inject Options ]");
    args::ValueFlag<std::string> _bed_targets(inject_opts, "FILE", "BED file over path space of the graph. Records will be converted into new paths in the output graph.",
                                              {'b', "bed-targets"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> nthreads(threading_opts, "N", "Number of threads to use for parallel operations.",
                                       {'t', "threads"});
    args::Group processing_info_opts(parser, "[ Processing Information ]");
    //args::Flag debug(processing_info_opts, "debug", "Print information about the process to stderr.", {'d', "debug"});
    args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi inject.", {'h', "help"});
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
            << "[odgi::inject] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
            << std::endl;
        return 1;
    }

    if (!og_out_file || args::get(og_out_file).empty()) {
        std::cerr << "[odgi::inject] error: please specify an output file to store the graph via -o=[FILE], --out=[FILE]." << std::endl;
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
                utils::handle_gfa_odgi_input(infile, "inject", args::get(progress), num_threads, graph);
            }
        }
    }

    ska::flat_hash_map<path_handle_t,
                       std::vector<std::pair<interval_t, std::string>>> path_intervals;
    std::vector<std::string> ordered_intervals;
    if (_bed_targets) {
        std::ifstream bed(args::get(_bed_targets).c_str());
        std::string line;
        while (std::getline(bed, line)) {
            if (!line.empty()) {
                auto vals = split(line, '\t');
                if (vals.size() < 3) {
                    std::cerr << "[odgi::inject]"
                              << "BED line does not have enough fields to define an interval"
                              << std::endl << line << std::endl;
                    return 1;
                }
                auto& path_name = vals[0];
                uint64_t start = std::stoul(vals[1]);
                uint64_t end = std::stoul(vals[2]);
                const std::string& name = vals[3];
                if (!graph.has_path(path_name)) {
                    //std::cerr << "[odgi::inject] warning: no path '" << path_name << "' in graph" << std::endl;
                } else {
                    auto path = graph.get_path_handle(path_name);
                    path_intervals[path].push_back(std::make_pair(interval_t(start, end), name));
                    ordered_intervals.push_back(name);
                }
            }
        }
        // ensure sorted input
        std::vector<std::vector<std::pair<interval_t, std::string>>*> v;
        for (auto& i : path_intervals) {
            v.push_back(&i.second);
        }
#pragma omp parallel for
        for (auto& ivals : v) {
            std::sort(ivals->begin(), ivals->end());
        }
    } else {
        std::cerr << "[odgi::inject] BED targets are required for injection of ranges" << std::endl;
        return 1;
    }

    omp_set_num_threads(num_threads);
    graph.set_number_of_threads(num_threads);

    algorithms::inject_ranges(graph, path_intervals, ordered_intervals, args::get(progress));

    const std::string outfile = args::get(og_out_file);
    if (outfile == "-") {
        graph.serialize(std::cout);
    } else {
        ofstream f(outfile.c_str());
        graph.serialize(f);
        f.close();
    }

    return 0;
}

static Subcommand odgi_inject("inject", "Inject BED annotations as paths.",
                              PIPELINE, 3, main_inject);


}
