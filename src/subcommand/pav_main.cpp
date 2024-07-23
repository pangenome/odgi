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
    --argc;

    args::ArgumentParser parser("Presence/absence variants (PAVs). It prints to stdout a TSV table with the 'PAV ratios'. "
                                "For a given path range 'PR' and path 'P', the 'PAV ratio' is the ratio between the sum of the lengths "
                                "of the nodes in 'PR' that are crossed by 'P' divided by the sum of the lengths"
                                " of all the nodes in 'PR'. Each node is considered only once.");
    args::Group mandatory_opts(parser, "[ MANDATORY ARGUMENTS ]");
    args::ValueFlag<std::string> og_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::ValueFlag<std::string> _path_bed_file(mandatory_opts, "FILE",
                                                "Find PAVs in the path range(s) specified in the given BED FILE.",
                                                {'b', "bed-file"});
    args::Group pav_opts(parser, "[ Pav Options ]");
    args::ValueFlag<std::string> _path_groups(pav_opts, "FILE", "Group paths as described in two-column FILE, with columns path.name and group.name.",
                                              {'p', "path-groups"});
    args::Flag _group_by_sample(pav_opts, "bool", "Following PanSN naming (sample#hap#ctg), group by sample (1st field).", {'S', "group-by-sample"});
    args::Flag _group_by_haplotype(pav_opts, "bool", "Following PanSN naming (sample#hap#ctg), group by haplotype (2nd field).", {'H', "group-by-haplotype"});
    args::ValueFlag<double> _binary_values(pav_opts, "THRESHOLD", "Print 1 if the PAV ratio is greater than or equal to the specified THRESHOLD, else 0.",
                                       {'B', "binary-values"});
    args::Flag _matrix_output(pav_opts, "bool", "Emit the PAV ratios in a matrix, with path ranges as rows and paths/groups as columns.",
                                           {'M', "matrix-output"});
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

    if (args::get(_binary_values) && (args::get(_binary_values) < 0 || args::get(_binary_values) > 1)) {
        std::cerr
            << "[odgi::pav] error: the PAV ratio threshold must be greater than 0 and lower than 1."
            << std::endl;
        return 1;
    }
    const double binary_threshold = args::get(_binary_values) ? args::get(_binary_values) : 0;
    const bool emit_binary_values = binary_threshold != 0;

    if (_path_groups + _group_by_sample + _group_by_haplotype > 1) {
        std::cerr
            << "[odgi::pav] error: select only one grouping option (-p/--path-groups, -S/--group-by-sample, or -H/--group-by-haplotype)."
            << std::endl;
        return 1;
    }

    // Read path groups
    const bool group_paths = _path_groups || _group_by_sample || _group_by_haplotype;
    ska::flat_hash_map<path_handle_t, std::string> path_2_group; // General, but potentially heavy solution
    std::map<std::string, uint64_t> group_2_index;               // Ordered map to keep group names' order
    if (group_paths) {
        if (_path_groups) {
            std::ifstream refs(args::get(_path_groups).c_str());
            std::string line;
            while (std::getline(refs, line)) {
                if (!line.empty()) {
                    auto vals = split(line, '\t');
                    if (vals.size() != 2) {
                        std::cerr << "[odgi::pav] line does not have a path.name and path.group value:"
                                  << std::endl << line << std::endl;
                        return 1;
                    }
                    auto& path_name = vals.front();
                    auto& group = vals.back();
                    if (!graph.has_path(path_name)) {
                        std::cerr << "[odgi::pav] no path '" << path_name << "'" << std::endl;
                        return 1;
                    }
                    path_2_group[graph.get_path_handle(path_name)] = group;
                    group_2_index[group] = 0;
                }
            }
            refs.close();

            if (group_2_index.empty()) {
                std::cerr
                        << "[odgi::pav] error: 0 path groups were read. Please specify at least one path group via -p/--path-groups."
                        << std::endl;
                return 1;
            }
        } else if (_group_by_sample) {
            graph.for_each_path_handle([&](const path_handle_t& p) {
                auto path_name = graph.get_path_name(p);
                // split and decide
                const auto vals = split(path_name, '#');
                const auto group = vals.front();
                path_2_group[p] = group;
                group_2_index[group] = 0;
            });
        } else if (_group_by_haplotype) {
            graph.for_each_path_handle([&](const path_handle_t& p) {
                auto path_name = graph.get_path_name(p);
                // split and decide
                const auto vals = split(path_name, '#');
                if (vals.size() == 1) { // Assume ctg
                    const auto group = vals.front();
                    path_2_group[p] = group;
                    group_2_index[group] = 0;
                } else if (vals.size() == 2) { // Assume sample#ctg
                    const auto group = vals.front();
                    path_2_group[p] = group;
                    group_2_index[group] = 0;
                } else if (vals.size() == 3) { // Assume sample#hap#ctg
                    const auto group = vals[0] + '#' + vals[1];
                    path_2_group[p] = group;
                    group_2_index[group] = 0;
                }
            });
        }

        uint64_t group_index = 0;
        for (auto& x : group_2_index) {
            x.second = group_index++;
        }
    }

    // Read target paths from BED
    std::vector<odgi::path_range_t> path_ranges;
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

    // Prepare the interval trees for querying target path ranges
    std::unique_ptr <odgi::algorithms::progress_meter::ProgressMeter> operation_progress;
    if (show_progress) {
        std::string banner = "[odgi::pav] preparing the interval trees for querying target path ranges:";
		operation_progress = std::make_unique<odgi::algorithms::progress_meter::ProgressMeter>(path_handles.size(), banner);
    }
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

        if (show_progress) {
            operation_progress->increment(1);
        }
    }
    if (show_progress) {
        operation_progress->finish();
    }
    const bool emit_matrix_else_table = args::get(_matrix_output);

    // Emit the PAV matrix
    std::cout << "chrom" << "\t"
              << "start" << "\t"
              << "end" << "\t"
              << "name";
    if (emit_matrix_else_table) {
        if (group_paths) {
            for (auto& x : group_2_index) {
                std::cout << "\t" << x.first;
            }
        } else {
            graph.for_each_path_handle([&](const path_handle_t path_handle) {
                std::cout << "\t" << graph.get_path_name(path_handle);
            });
        }
    } else {
        std::cout << "\t" << "group" << "\t" << "pav";
    }
    std::cout << std::endl;

    auto print_pav_table_row = [](
            const basic_ostream<char>& stream,
            graph_t& graph,
            const uint64_t len_unique_nodes_in_range,
            const std::vector<uint64_t>& len_unique_nodes_in_range_for_each_group,
            const std::string& group_name,
            const uint64_t group_rank,
            const odgi::path_range_t& path_range,
            const bool emit_binary_values,
            const double binary_threshold) {
        // Check if there were nodes in the range
        const double pav_ratio = len_unique_nodes_in_range == 0 ?
                                 0 : (double) len_unique_nodes_in_range_for_each_group[group_rank] / (double) len_unique_nodes_in_range;
        std::cout << std::setprecision(5)
                  << graph.get_path_name(path_range.begin.path) << "\t"
                  << path_range.begin.offset << "\t"
                  << path_range.end.offset << "\t"
                  << path_range.name << "\t"
                  << group_name << "\t"
                  << (emit_binary_values ? pav_ratio >= binary_threshold : pav_ratio) << "\n";
    };

    if (show_progress) {
        std::string banner = "[odgi::pav] emitting PAV results:";
		operation_progress = std::make_unique<odgi::algorithms::progress_meter::ProgressMeter>(path_ranges.size(), banner);
    }

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (uint64_t i = 0; i < path_ranges.size(); ++i) {
        auto &path_range = path_ranges[i];
        const uint64_t begin = path_range.begin.offset;
        const uint64_t end = path_range.end.offset;

        const uint64_t index = path_handle_2_index[path_range.begin.path];
        auto& tree = trees[index];

        std::vector<size_t> node_ids_info;
        tree.overlap(begin, end, node_ids_info); // retrieve overlaps

        uint64_t len_unique_nodes_in_range = 0;
        std::vector<uint64_t> len_unique_nodes_in_range_for_each_group(
                group_paths ? group_2_index.size() : graph.get_path_count(),
                0);

        // For each node in the range
        for (const auto& node_id_info : node_ids_info) {
            const auto& node_id = tree.data(node_id_info);
            const auto& handle = graph.get_handle(node_id);

            // Get paths that cross the node
            unordered_set<uint64_t> group_ranks_on_node_handle;
            graph.for_each_step_on_handle(handle, [&](const step_handle_t &source_step) {
                const auto& path_handle = graph.get_path_handle_of_step(source_step);
                // Check if the paths are grouped and there are paths that do not belong to any group
                if (!group_paths || path_2_group.find(path_handle) != path_2_group.end()) {
                    const uint64_t group_rank = group_paths ?
                            group_2_index[path_2_group[path_handle]] :
                            as_integer(path_handle) - 1;
                    group_ranks_on_node_handle.insert(group_rank);
                }
            });

            const uint64_t len_handle = graph.get_length(handle);
            for (const auto& group_rank: group_ranks_on_node_handle) {
                len_unique_nodes_in_range_for_each_group[group_rank] += len_handle;
            }

            len_unique_nodes_in_range += len_handle;
        }

#pragma omp critical (cout)
        {
            if (emit_matrix_else_table) {
                std::cout << std::setprecision(5)
                          << graph.get_path_name(path_range.begin.path) << "\t"
                          << path_range.begin.offset << "\t"
                          << path_range.end.offset << "\t"
                          << path_range.name;
                for (auto& x: len_unique_nodes_in_range_for_each_group) {
                    // Check if there were nodes in the range
                    const double pav_ratio = len_unique_nodes_in_range == 0 ?
                                             0 : (double) x / (double) len_unique_nodes_in_range;
                    std::cout << "\t" << (emit_binary_values ? pav_ratio >= binary_threshold : pav_ratio);
                }
                std::cout << std::endl;
            } else {
                if (group_paths) {
                    for (auto& x : group_2_index) {
                        const uint64_t group_rank = x.second;
                        print_pav_table_row(
                                std::cout,
                                graph,
                                len_unique_nodes_in_range,
                                len_unique_nodes_in_range_for_each_group,
                                x.first,
                                group_rank,
                                path_range,
                                emit_binary_values,
                                binary_threshold);
                    }
                } else {
                    graph.for_each_path_handle([&](const path_handle_t path_handle) {
                        const uint64_t group_rank = as_integer(path_handle) - 1;
                        print_pav_table_row(
                                std::cout,
                                graph,
                                len_unique_nodes_in_range,
                                len_unique_nodes_in_range_for_each_group,
                                graph.get_path_name(path_handle),
                                group_rank,
                                path_range,
                                emit_binary_values,
                                binary_threshold);
                    });
                }
            }
        }


        if (show_progress) {
            operation_progress->increment(1);
        }
    }
    if (show_progress) {
        operation_progress->finish();
    }
    return 0;
}

static Subcommand odgi_pav("pav", "Presence/absence variants (PAVs).",
                             PIPELINE, 3, main_pav);

}
