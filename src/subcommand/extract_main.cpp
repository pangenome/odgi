#include "subcommand.hpp"
#include "odgi.hpp"
#include "position.hpp"
#include "args.hxx"
#include "split.hpp"
#include <omp.h>
#include <regex>
#include "utils.hpp"
#include "atomic_bitvector.hpp"
#include "src/algorithms/subgraph/extract.hpp"
#include "path_length.hpp"

namespace odgi {

    using namespace odgi::subcommand;


    int main_extract(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        const std::string prog_name = "odgi extract";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser("Extract subgraphs or parts of a graph defined by query criteria.");
        args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
        args::ValueFlag<std::string> og_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
        args::Group graph_files_io_opts(parser, "[ Graph Files IO ]");
        args::ValueFlag<std::string> og_out_file(graph_files_io_opts, "FILE", "Store all subgraphs in this FILE. The file name usually ends with *.og*.",
                                                 {'o', "out"});
        args::Group extract_opts(parser, "[ Extract Options ]");
        args::ValueFlag<uint64_t> _max_dist_subpaths(extract_opts, "N",
                                                     "Maximum distance between subpaths allowed for merging them. "
                                                     "It reduces the fragmentation of unspecified paths in the input path ranges. "
                                                     "Set 0 to disable it [default: 300000].",
                                                     {'d', "max-distance-subpaths"});
        args::ValueFlag<uint64_t> _num_iterations(extract_opts, "N",
                                                  "Maximum number of iterations in attempting to merge close subpaths. "
                                                  "It stops early if during an iteration no subpaths were merged [default: 6].",
                                                  {'e', "max-merging-iterations"});
        args::Flag _split_subgraphs(extract_opts, "split_subgraphs",
                                    "Instead of writing the target subgraphs into a single graph, "
                                    "write one subgraph per given target to a separate file named path:start-end.og "
                                    "(0-based coordinates).", {'s', "split-subgraphs"});
        args::Flag _inverse(extract_opts, "inverse",
                               "Extract the parts of the graph that do not meet the query criteria.",
                               {'I', "inverse"});
        args::ValueFlag<uint64_t> _target_node(extract_opts, "ID", "A single node ID from which to begin our traversal.",
                                               {'n', "node"});
        args::ValueFlag<std::string> _node_list(extract_opts, "FILE", "A file with one node id per line. The node specified will be extracted from the input graph.", {'l', "node-list"});
        args::ValueFlag<uint64_t> _context_steps(extract_opts, "N",
                                                 "The number of steps (nodes) away from our initial subgraph that we should collect [default: 0 (disabled)]",
                                                 {'c', "context-steps"});
        args::ValueFlag<uint64_t> _context_bases(extract_opts, "N",
                                                 "The number of bases away from our initial subgraph that we should collect [default: 0 (disabled)]",
                                                 {'L', "context-bases"});
        args::ValueFlag<std::string> _path_range(extract_opts, "STRING",
                                                 "Find the node(s) in the specified path range TARGET=path[:pos1[-pos2]] "
                                                 "(0-based coordinates).", {'r', "path-range"});
        args::ValueFlag<std::string> _path_bed_file(extract_opts, "FILE",
                                                    "Find the node(s) in the path range(s) specified in the given BED FILE.",
                                                    {'b', "bed-file"});
        args::ValueFlag<std::string> _pangenomic_range(extract_opts, "STRING",
                                                       "Find the node(s) in the specified pangenomic range pos1-pos2 (0-based coordinates). "
                                                       "The nucleotide positions refer to the pangenome’s sequence (i.e., "
                                                       "the sequence obtained arranging all the graph’s node from left to right)."
                                                       , {'q', "pangenomic-range"});
        args::Flag _full_range(extract_opts, "full_range",
                               "Collects all nodes in the sorted order of the graph in the min and max positions touched by the given path ranges. "
                               "This ensures that all the paths of the subgraph are not split by node, but that the nodes are laced together again. "
                               "Comparable to **-R, --lace-paths=FILE**, but specifically for all paths in the resulting subgraph. "
                               "Be careful to use it with very complex graphs.",
                               {'E', "full-range"});
        args::ValueFlag<std::string> _path_names_file(extract_opts, "FILE",
                                                      "List of paths to keep in the extracted graph. The FILE must "
                                                      "contain one path name per line and a subset of all paths can be specified. "
                                                      "Paths specified in the input path ranges (with -r/--path-range and/or -b/--bed-file) "
                                                      "will be kept in any case.",
                                                      {'p', "paths-to-extract"});
        args::ValueFlag<std::string> _lace_paths_file(extract_opts, "FILE",
                                                       "List of paths to fully retain in the extracted graph. Must "
                                                       "contain one path name per line and a subset of all paths can be specified.",
                                                      {'R', "lace-paths"});
        args::Flag _optimize(extract_opts, "optimize", "Compact the node ID space in the extracted graph(s).",
                             {'O', "optimize"});
        args::Group threading_opts(parser, "[ Threading ]");
        args::ValueFlag<uint64_t> nthreads(threading_opts, "N", "Number of threads to use for parallel operations.",
                                           {'t', "threads"});
        args::Group processing_info_opts(parser, "[ Processing Information ]");
        args::Flag _show_progress(processing_info_opts, "progress",
                                  "Print information about the operations and the progress to stderr.",
                                  {'P', "progress"});
        args::Group program_info_opts(parser, "[ Program Information ]");
        args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi extract.", {'h', "help"});

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
                    << "[odgi::extract] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
                    << std::endl;
            return 1;
        }

        if ((_max_dist_subpaths && args::get(_max_dist_subpaths) == 0) && _num_iterations) {
            std::cerr << "[odgi::extract] error: specified -e/--max-merging-iterations without specifying -d/--max-distance-subpaths greater than 0." << std::endl;
            return 1;
        }

        if (_context_steps && _context_bases) {
            std::cerr << "[odgi::extract] error: please specify the expanding context either in steps (with -c/--context-steps) or "
                         "in bases (-L/--context-bases), not both." << std::endl;
            return 1;
        }

        if (args::get(_max_dist_subpaths) > 0 && _num_iterations && args::get(_num_iterations) == 0) {
            std::cerr << "[odgi::extract] error: -e/--max-merging-iterations has to be greater than 0." << std::endl;
            return 1;
        }

        const uint64_t max_dist_subpaths = _max_dist_subpaths && args::get(_max_dist_subpaths) >= 0 ? args::get(_max_dist_subpaths) : 300000;
        const uint64_t num_iterations = _num_iterations && args::get(_num_iterations) > 0 ? args::get(_num_iterations) : 6;

        if (_split_subgraphs) {
            if (og_out_file) {
                std::cerr << "[odgi::extract] error: please do not specify an output file (with -o/--out) when "
                             "one subgraph per given target is requested (with -s/--split-subgraphs)." << std::endl;
                return 1;
            }

            if (_target_node) {
                std::cerr << "[odgi::extract] error: please do not specify a single node (with -n/--node) when "
                             "one subgraph per given target is requested (with -s/--split-subgraphs)." << std::endl;
                return 1;
            }

            if (_node_list) {
                std::cerr << "[odgi::extract] error: please do not specify a node list (with -l/--list) when "
                             "one subgraph per given target is requested (with -s/--split-subgraphs)." << std::endl;
                return 1;
            }

            if (_inverse) {
                std::cerr << "[odgi::extract] error: please do not specify an inverse query (with -I/--inverse) when "
                             "one subgraph per given target is requested (with -s/--split-subgraphs)." << std::endl;
                return 1;
            }

            if (_pangenomic_range) {
                std::cerr << "[odgi::extract] error: please do not specify a pangenomic range (with -q/--pangenomic-range) when "
                             "one subgraph per given target is requested (with -s/--split-subgraphs)." << std::endl;
                return 1;
            }
        } else {
            if (!og_out_file) {
                std::cerr << "[odgi::extract] error: please specify an output file to where to store the subgraph via "
                             "-o=[FILE], --out=[FILE]." << std::endl;
                return 1;
            }
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
					utils::handle_gfa_odgi_input(infile, "extract", args::get(_show_progress), num_threads, graph);
                }
            }
        }

        const uint64_t shift = graph.min_node_id();
        if (graph.max_node_id() - shift >= graph.get_node_count()){
            std::cerr << "[odgi::extract] error: the node IDs are not compacted. Please run 'odgi sort' using -O, --optimize to optimize the graph." << std::endl;
            exit(1);
        }

        // Prepare all paths for parallelize the next step (actually, not all paths are always present in the subgraph)
        std::vector<path_handle_t> paths;
        if (args::get(_path_names_file).empty()) {
            paths.reserve(graph.get_path_count());
            graph.for_each_path_handle([&](const path_handle_t path) {
                paths.push_back(path);
            });
        } else {
            std::ifstream path_names_in(args::get(_path_names_file));

            uint64_t num_of_paths_in_file = 0;

            std::vector<bool> path_already_seen;
            path_already_seen.resize(graph.get_path_count(), false);

            std::string line;
            while (std::getline(path_names_in, line)) {
                if (!line.empty()) {
                    if (graph.has_path(line)) {
                        const path_handle_t path = graph.get_path_handle(line);
                        const uint64_t path_rank = as_integer(path) - 1;
                        if (!path_already_seen[path_rank]) {
                            path_already_seen[path_rank] = true;
                            paths.push_back(path);
                        } else {
                            std::cerr << "[odgi::extract] error: in the path list there are duplicated path names."
                                      << std::endl;
                            exit(1);
                        }
                    }

                    ++num_of_paths_in_file;
                }
            }

            path_names_in.close();

            std::cerr << "[odgi::extract] found " << paths.size() << "/" << num_of_paths_in_file
                      << " paths to consider." << std::endl;

            if (paths.empty()) {
                std::cerr << "[odgi::extract] error: no path to consider." << std::endl;
                exit(1);
            }
        }

        std::vector<path_handle_t> lace_paths;
        if (!args::get(_lace_paths_file).empty()) {
            ska::flat_hash_set<path_handle_t> lace_paths_set;
            std::ifstream path_names_in(args::get(_lace_paths_file));
            uint64_t num_of_paths_in_file = 0;
            std::string line;
            while (std::getline(path_names_in, line)) {
                if (!line.empty()
                    && graph.has_path(line)) {
                    auto path = graph.get_path_handle(line);
                    if (!lace_paths_set.count(path)) {
                        lace_paths.push_back(path);
                        lace_paths_set.insert(path);
                    }
                }
            }
            path_names_in.close();
            if (lace_paths.empty()) {
                std::cerr << "[odgi::extract] error: no path to fully retain." << std::endl;
                exit(1);
            }
        }

        std::vector<odgi::path_range_t> input_path_ranges;

        {
            // handle targets from BED
            if (_path_bed_file && !args::get(_path_bed_file).empty()) {
                std::ifstream bed_in(args::get(_path_bed_file));
                std::string line;
                while (std::getline(bed_in, line)) {
                    add_bed_range(input_path_ranges, graph, line);
                }
            }

            // handle targets from command line
            if (_path_range) {
                Region region;
                parse_region(args::get(_path_range), region);

                // no coordinates given, we do whole thing (0,-1)
                if (region.start < 0 || region.end < 0) {
                    add_bed_range(input_path_ranges, graph, region.seq);
                } else {
                    add_bed_range(input_path_ranges, graph, region.seq + "\t" + std::to_string(region.start) + "\t" + std::to_string(region.end));
                }
            }

            // Check duplicates
            std::vector<path_range_t> copy_ranges = input_path_ranges; // Create a copy of the vector to avoid sorting the original one

            auto compare_path_range = [](const path_range_t& a, const path_range_t& b) -> bool {
                if (a.begin.path != b.begin.path) return a.begin.path < b.begin.path;
                if (a.begin.offset != b.begin.offset) return a.begin.offset < b.begin.offset;
                if (a.end.path != b.end.path) return a.end.path < b.end.path;
                return a.end.offset < b.end.offset;
            }; // Lambda function to compare two path_range_t objects

            std::sort(copy_ranges.begin(), copy_ranges.end(), compare_path_range); // Sort the copied vector using the lambda function

            for (size_t i = 1; i < copy_ranges.size(); i++) {
                if (!compare_path_range(copy_ranges[i-1], copy_ranges[i])) {
                    std::cerr << "[odgi::extract] error: " << graph.get_path_name(copy_ranges[i].begin.path) << ":" << copy_ranges[i].begin.offset << "-" << copy_ranges[i].end.offset << " is a duplicated path range" << std::endl;
                    return 1;
                }
            }
        }

        std::vector<std::pair<uint64_t, uint64_t>> input_pangenomic_ranges;
        if (_pangenomic_range && !args::get(_pangenomic_range).empty()) {
            const std::string pangenomic_range_str = args::get(_pangenomic_range);

            const std::regex regex("-");
            const std::vector<std::string> splitted(
                    std::sregex_token_iterator(pangenomic_range_str.begin(), pangenomic_range_str.end(), regex, -1),
                    std::sregex_token_iterator()
                    );

            if (splitted.size() != 2) {
                std::cerr
                        << "[odgi::extract] error: please specify a valid pangenomic range: start-end."
                        << std::endl;
                return 1;
            }

            if (!utils::is_number(splitted[0]) || !utils::is_number(splitted[1]) || stoull(splitted[0]) > stoull(splitted[1])) {
                std::cerr
                << "[odgi::extract] error: please specify valid numbers for the pangenomic range."
                << std::endl;
                return 1;
            }

            input_pangenomic_ranges.push_back({stoull(splitted[0]), stoull(splitted[1])});
        }

        if (_split_subgraphs && input_path_ranges.empty()) {
            std::cerr << "[odgi::extract] error: please specify at least one target when "
                         "one subgraph per given target is requested (with -s/--split-subgraphs)." << std::endl;
            exit(1);
        }

        // Path ranges inversion
        std::vector<odgi::path_range_t>* path_ranges;
        if (_inverse && !input_path_ranges.empty()) {
            auto path_length_table = algorithms::get_path_length(graph);

            // On average, each interval (in the middle) would generate two intervals (on the sides) when inverted
            path_ranges = new std::vector<odgi::path_range_t>();
            path_ranges->reserve(input_path_ranges.size() * 2);

            for (auto &path_range : input_path_ranges) {
                const path_handle_t path_handle = path_range.begin.path;
                const uint64_t path_length = path_length_table[path_handle];

                if (path_range.begin.offset > 0) {
                    path_ranges->push_back({
                        { path_handle, 0, false },
                        { path_handle, path_range.begin.offset, false },
                        path_range.is_rev,
                        path_range.name == "." ? "." : path_range.name + "-1",
                        ""
                    });
                }

                if (path_range.end.offset < path_length) {
                    path_ranges->push_back({
                        { path_handle, path_range.end.offset, false },
                        { path_handle, path_length, false },
                        path_range.is_rev,
                        path_range.name == "." ? "." : path_range.name + "-2",
                        ""
                    });
                }
            }

            input_path_ranges.clear();
        } else {
            path_ranges = &input_path_ranges;
        }

        // Pangenomic range inversion
        std::vector<std::pair<uint64_t, uint64_t>>* pangenomic_ranges;
        if (_inverse && _pangenomic_range) {
            uint64_t pangenome_len = 0;
            graph.for_each_handle([&](const handle_t &h) {
                pangenome_len += graph.get_length(h);
            });

            pangenomic_ranges = new std::vector<std::pair<uint64_t, uint64_t>>();

            if (input_pangenomic_ranges[0].first > 0) {
                pangenomic_ranges->push_back({0, input_pangenomic_ranges[0].first});
            }
            if (input_pangenomic_ranges[0].second < pangenome_len) {
                pangenomic_ranges->push_back({input_pangenomic_ranges[0].second, pangenome_len});
            }
        } else {
            pangenomic_ranges = &input_pangenomic_ranges;
        }


        const bool show_progress = args::get(_show_progress);
        const bool optimize = args::get(_optimize);
        const uint64_t context_steps = _context_steps ? args::get(_context_steps) : 0;
        const uint64_t context_bases = _context_bases ? args::get(_context_bases) : 0;

        omp_set_num_threads((int) num_threads);

        auto prep_graph = [&shift](
                             graph_t &source, std::vector<path_handle_t>* source_paths,
                             const std::vector<path_handle_t>& lace_paths, graph_t &subgraph,
                             std::vector<odgi::path_range_t> path_ranges, std::vector<std::pair<uint64_t, uint64_t>> pangenomic_ranges,
                             const uint64_t context_steps, const uint64_t context_bases, const bool full_range, const bool inverse,
                             const uint64_t max_dist_subpaths, const uint64_t num_iterations,
                             const uint64_t num_threads, const bool show_progress, const bool optimize) {
            if (context_steps > 0 || context_bases > 0) {
                if (show_progress) {
                    std::cerr << "[odgi::extract] expansion and adding connecting edges" << std::endl;
                }

                if (context_steps > 0) {
                    algorithms::expand_subgraph_by_steps(source, subgraph, context_steps, false);
                } else {
                    algorithms::expand_subgraph_by_length(source, subgraph, context_bases, false);
                }
            }

            // Check if there are nodes in the subgraph, to avoid extracting the whole graph
            // when nodes with functionality other than pangenomic paths/ranges are not specified
            if (inverse && subgraph.get_node_count() > 0) {
                unordered_set<nid_t> node_ids_to_ignore;
                subgraph.for_each_handle([&](const handle_t &h) {
                    node_ids_to_ignore.insert(subgraph.get_id(h));
                });

                std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress;
                if (show_progress) {
                    progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                            source.get_node_count() - node_ids_to_ignore.size(), "[odgi::extract] inverting the query criteria");
                }

                subgraph.clear();

                source.for_each_handle([&](const handle_t &h) {
                    nid_t id = source.get_id(h);
                    if (node_ids_to_ignore.count(id) <= 0) {
                        subgraph.create_handle(source.get_sequence(source.get_handle(id)), id);

                        if (show_progress) {
                            progress->increment(1);
                        }
                    }
                });

                if (show_progress) {
                    progress->finish();
                }
            }

            // Collect handles in path/pangenomic ranges (it is assumed they were already inverted outside, if needed)
            {
                std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress;
                if (show_progress) {
                    progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                            path_ranges.size(), "[odgi::extract] extracting path ranges");
                }

                atomicbitvector::atomic_bv_t keep_bv(source.get_node_count()+1);

#pragma omp parallel for schedule(dynamic,1)
                for (auto &path_range : path_ranges) {
                    const path_handle_t path_handle = path_range.begin.path;
                    if (show_progress) {
                        progress->increment(1);
                    }

                    // The extraction does not cut nodes, so the input path ranges have to be
                    // extended if their ranges (start, end) fall in the middle of the nodes.
                    bool first = true;
                    uint64_t new_start = 0;
                    uint64_t new_end = 0;
                    
                    const uint64_t start = path_range.begin.offset;
                    const uint64_t end = path_range.end.offset;

                    uint64_t walked = 0;
                    const auto path_end = source.path_end(path_handle);
                    for (step_handle_t cur_step = source.path_begin(path_handle);
                        cur_step != path_end && walked < end; cur_step = source.get_next_step(cur_step)) {
                        const handle_t cur_handle = source.get_handle_of_step(cur_step);
                        walked += source.get_length(cur_handle);
                        if (walked > start) {
                            keep_bv.set(source.get_id(cur_handle) - shift);

                            if (first) {
                                first = false;
                                new_start = walked - source.get_length(cur_handle);
                            }
                        }
                    }
                    new_end = walked;

                    // Extend path range to entirely include the first and the last node of the range.
                    // Thi is important to path names with the correct path ranges.
                    path_range.begin.offset = new_start;
                    path_range.end.offset = new_end;
                }
                if (!pangenomic_ranges.empty()) {
                    uint64_t pos = 0;
                    source.for_each_handle([&](const handle_t &h) {
                        const uint64_t hl = source.get_length(h);

                        for (auto &pan_range : pangenomic_ranges) {
                            if (pos + hl >= pan_range.first && pos <= pan_range.second ) {
                                keep_bv.set(source.get_id(h) - shift);
                                break;
                            }
                        }

                        pos += hl;
                    });
                }
                for (auto id_shifted : keep_bv) {
                    const handle_t h = source.get_handle(id_shifted + shift);
                    subgraph.create_handle(source.get_sequence(h),
                                           id_shifted + shift);
                }

                if (show_progress) {
                    progress->finish();
                }
            }

            // Check if there are nodes in the subgraph, to avoid min_node_id == max_node_id == 0
            if (full_range && subgraph.get_node_count() > 0) {
                // Take the start and end node of this and fill things in
                algorithms::extract_id_range(source, subgraph.min_node_id(), subgraph.max_node_id() , subgraph, show_progress
                                                                                                                ? "[odgi::extract] collecting all nodes in the path range"
                                                                                                                : "");
            }

            // These paths are treated differently: only the specified ranges are included in the extracted graph
            std::vector<path_handle_t> source_paths_from_path_ranges;
            for (auto &path_range : path_ranges) {
                source_paths_from_path_ranges.push_back(path_range.begin.path);
            }

            // `max_dist_subpaths` and `add_subpaths_to_subgraph` have to work with the paths not specified in the
            // input path ranges, preventing their possible fragmentation.
            std::sort(source_paths_from_path_ranges.begin(), source_paths_from_path_ranges.end());
            source_paths->erase(std::remove_if(source_paths->begin(), source_paths->end(), [&](const auto&x) {
                return std::binary_search(source_paths_from_path_ranges.begin(), source_paths_from_path_ranges.end(), x);
            }), source_paths->end());

            // We don't cut nodes for the extraction, so close path intervals can generate identical subpaths.
            // To avoid duplicated subpaths in the final subgraph, we remove duplicated path ranges.
            {
                std::set<odgi::path_range_t, odgi::path_range_comparator> unique_path_ranges;

                for (const auto& path_range : path_ranges) {
                    unique_path_ranges.insert(path_range);
                }

                path_ranges.assign(unique_path_ranges.begin(), unique_path_ranges.end());
            }

            if (max_dist_subpaths > 0) {
                // Iterate multiple times to merge subpaths which became mergeable during the first iteration where new nodes were added
                for (uint8_t i = 0; i < num_iterations; ++i) {
                    std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress;
                    if (show_progress) {
                        progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                                source.get_path_count(), "[odgi::extract] merge subpaths closer than " + std::to_string(max_dist_subpaths) + " bps - iteration " +
                                                         std::to_string(i + 1) + " (max " + std::to_string(num_iterations) + ")");
                    }

                    // The last step is not included
                    std::vector<std::pair<step_handle_t, step_handle_t>> short_missing_subpaths;

                    // Search not included subpaths (in parallel)
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
                    for (uint64_t path_rank = 0; path_rank < source_paths->size(); ++path_rank) {
                        auto &source_path_handle = (*source_paths)[path_rank];
                        uint64_t walked = 0;

                        // check if the nodes are in the output subgraph
                        bool in_match = true;
                        step_handle_t start_step = source.path_begin(source_path_handle);

                        // We don't have to consider the not included subpath at the beginning (if any)
                        bool ignore_subpath = !subgraph.has_node(source.get_id(source.get_handle_of_step(start_step)));

                        uint64_t start_nt, end_nt;
                        // get each range that isn't included
                        source.for_each_step_in_path(source_path_handle, [&](const step_handle_t& step) {
                            const handle_t source_handle = source.get_handle_of_step(step);
                            const uint64_t source_length = source.get_length(source_handle);

                            if (!subgraph.has_node(source.get_id(source_handle))) {
                                if (in_match) {
                                    in_match = false;
                                    start_step = step;
                                    start_nt = walked;
                                }

                                end_nt = walked + source_length;
                            } else {
                                if (!in_match) {
                                    if (!ignore_subpath) {
                                        if ((end_nt - start_nt) <= max_dist_subpaths) {
#pragma omp critical (short_missing_subpaths)
                                            short_missing_subpaths.push_back(std::make_pair(start_step, step));
                                        }
                                    }

                                    ignore_subpath = false;
                                }
                                in_match = true;
                            }

                            walked += source_length;
                        });
                        // Ignore last not included subpath (if any)

                        if (show_progress) {
                            progress->increment(1);
                        }
                    }

                    if (show_progress) {
                        progress->finish();
                    }

                    if (short_missing_subpaths.empty()) {
                        break; // Nothing mergeable, do not waste time in further iterations
                    }

                    // Restore short subpaths by adding the associated handles
                    for (auto& range : short_missing_subpaths) {
                        if (range.first != range.second) {
                            for (step_handle_t step = range.first; step != range.second; step = source.get_next_step(step)) {
                                handle_t h = source.get_handle_of_step(step);
                                const uint64_t id = source.get_id(h);
                                // To avoid adding multiple times the same node
                                if(!subgraph.has_node(id)) {
                                    if (source.get_is_reverse(h)) {
                                        h = source.flip(h); // All handles are added in forward in the subgraph
                                    }
                                    subgraph.create_handle(source.get_sequence(h), id);
                                }
                            }
                        }
                    }
                }
            }

            // ----------------------------------------------------------------------------------
            // Insert the subpaths corresponding to the path ranges (if any)
            // Create subpaths
            std::vector<path_handle_t> subpaths_from_path_ranges;
            subpaths_from_path_ranges.reserve(path_ranges.size());

            for (auto &path_range : path_ranges) {
                const std::string path_name = source.get_path_name(path_range.begin.path);

                subpaths_from_path_ranges.push_back(
                        // The function assumes that every path is new and unique
                        odgi::algorithms::create_subpath(
                            subgraph,
                            odgi::algorithms::make_path_name(path_name, path_range.begin.offset, path_range.end.offset),
                            source.get_is_circular(path_range.begin.path)
                    )
                );
            }

            // Fill subpaths in parallel
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
            for (uint64_t i = 0; i < subpaths_from_path_ranges.size(); ++i) {
                auto path_range = path_ranges[i];
                const path_handle_t path_handle = path_range.begin.path;
                const path_handle_t subpath_handle = subpaths_from_path_ranges[i];

                algorithms::for_handle_in_path_range(
                        source, path_handle, path_range.begin.offset, path_range.end.offset,
                        [&](const handle_t& handle) {
                            subgraph.append_step(
                                    subpath_handle,
                                    subgraph.get_handle(source.get_id(handle),
                                                        source.get_is_reverse(handle))
                            );
                        });
            }
            // ----------------------------------------------------------------------------------

            // rewrite lace paths so that skipped regions are represented as new nodes that we then add to our subgraph
            if (!lace_paths.empty()) {
                if (show_progress) {
                    std::cerr << "[odgi::extract] adding " << lace_paths.size() << " lace paths" << std::endl;
                }

                algorithms::embed_lace_paths(source, subgraph, lace_paths);
            }

            // Connect the collected handles
            algorithms::add_connecting_edges_to_subgraph(source, subgraph, show_progress
                                                                           ? "[odgi::extract] adding connecting edges"
                                                                           : "");

            // Add subpaths covering the collected handles
            algorithms::add_subpaths_to_subgraph(source, *source_paths, subgraph, num_threads,
                                                 show_progress ? "[odgi::extract] adding subpaths" : "");

            std::vector<path_handle_t> subpaths;
            subpaths.reserve(subgraph.get_path_count());
            subgraph.for_each_path_handle([&](const path_handle_t& path) {
                subpaths.push_back(path);
            });

            std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_checking;
            if (show_progress) {
                progress_checking = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                        subpaths.size(), "[odgi::extract] checking missing edges and empty subpaths");
            }

            ska::flat_hash_set<std::pair<handle_t, handle_t>> edges_to_create;

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
            for (auto path: subpaths) {
                handle_t last;
                const step_handle_t begin_step = subgraph.path_begin(path);
                subgraph.for_each_step_in_path(path, [&](const step_handle_t &step) {
                    handle_t h = subgraph.get_handle_of_step(step);
                    if (step != begin_step && !subgraph.has_edge(last, h)) {
#pragma omp critical (edges_to_create)
                        edges_to_create.insert({last, h});
                    }
                    last = h;
                });

                if (show_progress) {
                    progress_checking->increment(1);
                }
            }

            if (show_progress) {
                progress_checking->finish();
            }

            // remove empty subpaths
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
            for (auto path: subpaths) {
                if (subgraph.is_empty(path)) {
#pragma omp critical (subgraph)
                    subgraph.destroy_path(path);
                }
            }

            if (show_progress && subgraph.get_path_count() < subpaths.size()) {
                std::cerr << "[odgi::extract] removed " << (subpaths.size() - subgraph.get_path_count()) << " empty subpath(s)." << std::endl;
            }

            subpaths.clear();

            // add missing edges
            for (auto edge: edges_to_create) {
                subgraph.create_edge(edge.first, edge.second);
            }

            if (show_progress && edges_to_create.size() > 0) {
                std::cerr << "[odgi::extract] fixed " << edges_to_create.size() << " edge(s)" << std::endl;
            }

            // This should not be necessary, if the extraction works correctly
            // subgraph.remove_orphan_edges();

            if (optimize) {
                subgraph.optimize();
            }
        };

        auto check_and_create_handle = [&](const graph_t &source, graph_t &subgraph, const nid_t node_id) {
            if (graph.has_node(node_id)) {
                if (!subgraph.has_node(node_id)){
                    const handle_t cur_handle = graph.get_handle(node_id);
                    subgraph.create_handle(
                            source.get_sequence(source.get_is_reverse(cur_handle) ? source.flip(cur_handle) : cur_handle),
                            node_id);
                }
            } else {
                std::cerr << "[odgi::extract] warning, cannot find node " << node_id << std::endl;
            }
        };

        if (_split_subgraphs) {
            for (auto &path_range : *path_ranges) {
                graph_t subgraph;

                const path_handle_t path_handle = path_range.begin.path;

                if (show_progress) {
                    std::cerr << "[odgi::extract] extracting path range " << graph.get_path_name(path_range.begin.path) << ":" << path_range.begin.offset
                              << "-"
                              << path_range.end.offset << std::endl;
                }

                prep_graph(
                    graph, &paths,
                    lace_paths, subgraph,
                    {path_range}, *pangenomic_ranges,
                    context_steps, context_bases, _full_range, false,
                    max_dist_subpaths, num_iterations,
                    num_threads, show_progress, optimize);

                const string filename = graph.get_path_name(path_range.begin.path) + ":" + to_string(path_range.begin.offset) + "-" + to_string(path_range.end.offset) + ".og";

                if (show_progress) {
                    std::cerr << "[odgi::extract] writing " << filename << std::endl;
                }

                ofstream f(filename);
                subgraph.serialize(f);
                f.close();
            }
        } else {
            graph_t subgraph;

            if (args::get(_target_node)) {
                check_and_create_handle(graph, subgraph, args::get(_target_node));
            }

            if (!args::get(_node_list).empty()) {
                ifstream nodes(args::get(_node_list));
                std::string next;
                while (std::getline(nodes, next)) {
                    check_and_create_handle(graph, subgraph, std::stol(next));
                }
            }

            prep_graph(
                graph, &paths,
                lace_paths, subgraph,
                *path_ranges, *pangenomic_ranges,
                context_steps, context_bases, _full_range, _inverse,
                max_dist_subpaths, num_iterations,
                num_threads, show_progress, optimize);

            {
                const std::string outfile = args::get(og_out_file);
                if (!outfile.empty()) {
                    if (outfile == "-") {
                        subgraph.serialize(std::cout);
                    } else {
                        ofstream f(outfile.c_str());
                        subgraph.serialize(f);
                        f.close();
                    }
                }
            }
        }

        return 0;
    }

    static Subcommand odgi_extract("extract", "Extract subgraphs or parts of a graph defined by query criteria.",
                                   PIPELINE, 3, main_extract);

}
