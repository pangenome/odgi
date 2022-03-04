#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include <position.hpp>
#include <subgraph/extract.hpp>
#include "utils.hpp"
#include "split.hpp"
#include "subgraph/region.hpp"
#include "IITree.h"

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
    args::ArgumentParser parser("Presence/absence variants (PAVs). It prints to stdout a matrix with the PAVs ratios.");
    args::Group mandatory_opts(parser, "[ MANDATORY ARGUMENTS ]");
    args::ValueFlag<std::string> og_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::Group pav_opts(parser, "[ Pav Options ]");
    args::ValueFlag<std::string> _path_bed_file(pav_opts, "FILE",
                                                "Find PAVs in the path range(s) specified in the given BED FILE.",
                                                {'b', "bed-file"});
    args::ValueFlag<double> _binary_matrix(pav_opts, "THRESHOLD", "Emit a binary matrix, with 1 if the PAV ratio is greater than or equal to the specified THRESHOLD, else 0.",
                                       {'B', "binary-matrix"});
//    args::ValueFlag<std::string> _path_groups(pav_opts, "FILE", "Group paths as described in two-column FILE, with columns path.name and group.name.",
//                                              {'p', "path-groups"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> nthreads(threading_opts, "N", "Number of threads to use for parallel operations.",
                                       {'t', "threads"});
    args::Group processing_info_opts(parser, "[ Processing Information ]");
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

    if (args::get(_binary_matrix) && (args::get(_binary_matrix) < 0 || args::get(_binary_matrix) > 1)) {
        std::cerr
            << "[odgi::pav] error: the PAV ratio must be greather than 0 and lower than 1."
            << std::endl;
        return 1;
    }

    const double binary_threshold = args::get(_binary_matrix) ? args::get(_binary_matrix) : 0;
    const bool emit_binary_matrix = binary_threshold != 0;

    std::vector<odgi::path_range_t> path_ranges;

    // Read target paths from BED
    if (_path_bed_file && !args::get(_path_bed_file).empty()) {
        std::ifstream bed_in(args::get(_path_bed_file));
        std::string line;
        while (std::getline(bed_in, line)) {
            add_bed_range(path_ranges, graph, line);
        }
    }

    if (path_ranges.empty()) {
        std::cerr
            << "[odgi::pav] error: please specify path ranges via -b/--bed-file."
            << std::endl;
        return 1;
    }

    std::vector<unordered_set<handle_t>> path_range_node_handles;
    path_range_node_handles.resize(path_ranges.size());

    // Check the min/max coordinates for each target path
    ska::flat_hash_map<path_handle_t, std::pair<uint64_t, uint64_t>> path_name_2_min_max;
    for (uint64_t i = 0; i < path_ranges.size(); ++i) {
        auto &path_range = path_ranges[i];
        path_name_2_min_max[path_range.begin.path] = { std::numeric_limits<uint64_t>::max(), 0 };
    }
    for (uint64_t i = 0; i < path_ranges.size(); ++i) {
        auto &path_range = path_ranges[i];
        const uint64_t begin = path_range.begin.offset;
        const uint64_t end = path_range.end.offset;

        if (begin < path_name_2_min_max[path_range.begin.path].first) {
            path_name_2_min_max[path_range.begin.path].first = begin;
        }
        if (end > path_name_2_min_max[path_range.begin.path].second) {
            path_name_2_min_max[path_range.begin.path].second = end;
        }
    }

    // Prepare target path handles to run in parallel on them
    ska::flat_hash_map<path_handle_t, uint64_t> path_handle_2_index;
    std::vector<path_handle_t> path_handles;
    path_handles.reserve(path_name_2_min_max.size());
    uint64_t i = 0;
    for (auto& p2mm : path_name_2_min_max) {
        path_handles.push_back(p2mm.first);
        path_handle_2_index[p2mm.first] = i++;
    }

    // Prepare the interval trees to query target path ranges
    std::vector<IITree<uint64_t, uint64_t>> trees;
    trees.resize(path_handles.size());
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (uint64_t i = 0; i < path_handles.size(); ++i) {
        const auto& path_handle = path_handles[i];
        const uint64_t min = path_name_2_min_max[path_handle].first;
        const uint64_t max = path_name_2_min_max[path_handle].second;

        const uint64_t index = path_handle_2_index[path_handle];
        auto& tree = trees[index];

        uint64_t walked = 0;
        const auto path_end = graph.path_end(path_handle);
        for (step_handle_t cur_step = graph.path_begin(path_handle);
        cur_step != path_end && walked < max; cur_step = graph.get_next_step(cur_step)) {
            const handle_t cur_handle = graph.get_handle_of_step(cur_step);
            const uint64_t len_cur_handle = graph.get_length(cur_handle);
            walked += len_cur_handle;
            if (walked > min) {
                tree.add(walked - len_cur_handle, walked, graph.get_id(cur_handle));
            }
        }

        tree.index(); // index
    }

    // Emit the PAV matrix
    //todo we are not managing the strandness... should we?
    std::cout << "chrom" << "\t"
    << "start" << "\t"
    << "end";
    graph.for_each_path_handle([&](const path_handle_t path_handle) {
        std::cout << "\t" << graph.get_path_name(path_handle);
    });
    std::cout << std::endl;

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (uint64_t i = 0; i < path_ranges.size(); ++i) {
        auto &path_range = path_ranges[i];
        const uint64_t begin = path_range.begin.offset;
        const uint64_t end = path_range.end.offset;

        const uint64_t index = path_handle_2_index[path_range.begin.path];
        auto& tree = trees[index];

        std::vector<uint64_t> node_ids_info;
        tree.overlap(begin, end, node_ids_info); // retrieve overlaps
        if (!node_ids_info.empty()) {
            uint64_t len_unique_nodes_in_range = 0;
            std::vector<uint64_t> len_unique_nodes_in_range_each_path(graph.get_path_count(), 0);

            // For each node in the range
            for (const auto& node_id_info : node_ids_info) {
                const auto& node_id = tree.data(node_id_info);
                const auto& handle = graph.get_handle(node_id);

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
                    const double pav_ratio = (double) x / (double) len_unique_nodes_in_range;
                    std::cout << "\t" << (emit_binary_matrix ? pav_ratio >= binary_threshold : pav_ratio);
                }
                std::cout << std::endl;
            }
        }
    }

    return 0;
}

static Subcommand odgi_pav("pav", "Presence/absence variants (PAVs).",
                             PIPELINE, 3, main_pav);

}
