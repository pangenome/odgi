#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include <position.hpp>
#include <subgraph/extract.hpp>
#include "utils.hpp"
#include "split.hpp"
#include "subgraph/region.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_pav(int argc, char **argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc - 1; ++i) {
        argv[i] = argv[i + 1];
    }
    const std::string prog_name = "odgi pav";
    argv[0] = (char *) prog_name.c_str();
    --argc
;
    args::ArgumentParser parser("Presence/absence variants (PAVs).");
    args::Group mandatory_opts(parser, "[ MANDATORY ARGUMENTS ]");
    args::ValueFlag<std::string> og_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::Group pav_opts(parser, "[ Pav Options ]");
    args::ValueFlag<std::string> _path_bed_file(pav_opts, "FILE",
                                                "Find PAVs in the path range(s) specified in the given BED FILE.",
                                                {'b', "bed-file"});
    args::ValueFlag<std::string> _path_groups(pav_opts, "FILE", "Group paths as described in two-column FILE, with columns path.name and group.name.",
                                              {'p', "path-groups"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> nthreads(threading_opts, "N", "Number of threads to use for parallel operations.",
                                       {'t', "threads"});
    args::Group processing_info_opts(parser, "[ Processing Information ]");
    args::Flag debug(processing_info_opts, "debug", "Print information about the process to stderr.", {'d', "debug"});
    args::Flag _progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi pav.", {'h', "help"});

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
        << "[odgi::pav] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
            << std::endl;
        return 1;
    }

    const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;
    omp_set_num_threads(num_threads);

    const bool show_progress = args::get(_progress);

    graph_t graph;
    assert(argc > 0);
    {
        const std::string infile = args::get(og_in_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
                utils::handle_gfa_odgi_input(infile, "pav", show_progress, num_threads, graph);
            }
        }
    }
    graph.set_number_of_threads(num_threads);

    std::vector<odgi::path_range_t> path_ranges;

    // handle targets from BED
    if (_path_bed_file && !args::get(_path_bed_file).empty()) {
        std::ifstream bed_in(args::get(_path_bed_file));
        std::string line;
        while (std::getline(bed_in, line)) {
            add_bed_range(path_ranges, graph, line);
        }
    }

    std::vector<unordered_set<handle_t>> path_range_node_handles;
    path_range_node_handles.resize(path_ranges.size());

    if (!path_ranges.empty()) {
        std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress;
        if (show_progress) {
            progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                    path_ranges.size(), "[odgi::pav] collect node handles for the path range");
        }
        #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
        for (uint64_t i = 0; i < path_ranges.size(); ++i) {
            auto &path_range = path_ranges[i];

            algorithms::for_handle_in_path_range(
                    graph, path_range.begin.path, path_range.begin.offset, path_range.end.offset,
                    [&](const handle_t& cur_handle) {
                        path_range_node_handles[i].insert(cur_handle);
                    });
            if (show_progress) {
                progress->increment(1);
            }
        }
        if (show_progress) {
            progress->finish();
        }
    }

    // Prepare all paths for parallelize the next step
    // todo to support subsets?
//    std::vector<path_handle_t> paths;
//    paths.reserve(graph.get_path_count());
//    graph.for_each_path_handle([&](const path_handle_t path) {
//        paths.push_back(path);
//    });

    //todo we are not managing the strandness... should we?
    //todo we are still managing entire nodes, not chunks
    std::cout << "chrom" << "\t"
              << "start" << "\t"
              << "end";
    graph.for_each_path_handle([&](const path_handle_t path_handle) {
        std::cout << graph.get_path_name(path_handle) << "\t";
    });
    std::cout << std::endl;

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (uint64_t i = 0; i < path_ranges.size(); ++i) {
        auto &path_range = path_ranges[i];

        uint64_t len_unique_nodes_in_range = 0;
        std::vector<uint64_t> len_unique_nodes_in_range_each_path(graph.get_path_count(), 0);

        // For each node in the range
        for (const auto& handle : path_range_node_handles[i]) {
            // Get paths that cross the node
            unordered_set<uint64_t> path_handles_on_handle;
            graph.for_each_step_on_handle(handle, [&](const step_handle_t &source_step) {
                const auto& path_handle = graph.get_path_handle_of_step(source_step);
                const uint64_t path_rank = as_integer(path_handle) - 1;
                path_handles_on_handle.insert(path_rank);
            });

            const uint64_t len_handle = graph.get_length(handle);
            for (const auto& path_rank: path_handles_on_handle) {
                len_unique_nodes_in_range_each_path[path_rank] += len_handle;
            }

            len_unique_nodes_in_range += len_handle;
        }

#pragma omp critical (cout)
        {
            std::cout << std::setprecision(5)
                      << graph.get_path_name(path_range.begin.path) << "\t"
                      << path_range.begin.offset << "\t"
                      << path_range.end.offset;
            for (auto& x: len_unique_nodes_in_range_each_path) {
                std::cout << "\t" << (double) x / (double) len_unique_nodes_in_range;
            }
            std::cout << std::endl;
        }

//        if (path_handle == path_range.begin.path) {
//            continue; // Skip itself
//        }
//#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
//        for (uint64_t path_rank = 0; path_rank < paths.size(); ++path_rank) {
    }

    return 0;
}

static Subcommand odgi_pav("pav", "Presence/absence variants (PAVs).",
                             PIPELINE, 3, main_pav);


}
