#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/prune.hpp"
#include "algorithms/depth.hpp"
#include "algorithms/remove_high_degree.hpp"
#include "algorithms/cut_tips.hpp"
#include "algorithms/remove_isolated.hpp"
#include "algorithms/expand_context.hpp"
#include "utils.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_prune(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    const std::string prog_name = "odgi prune";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("Remove parts of the graph.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1 or GFAz (compressed GFA), but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::ValueFlag<std::string> dg_out_file(mandatory_opts, "FILE", "Write the pruned graph in ODGI format to *FILE*. A file ending with *.og* is recommended.", {'o', "out"});
    args::Group kmer_opts(parser, "[ Kmer Options ]");
    args::ValueFlag<uint64_t> kmer_length(kmer_opts, "K", "The length of the kmers to consider.", {'k', "kmer-length"});
    args::ValueFlag<uint64_t> max_furcations(kmer_opts, "N", "Break at edges that would induce *N* many furcations in a kmer.", {'e', "max-furcations"});
    args::Group node_options(parser, "[ Node Options ]");
    args::ValueFlag<uint64_t> max_degree(node_options, "N", "Remove nodes that have a higher node degree than *N*.", {'d', "max-degree"});
    args::ValueFlag<uint64_t> min_depth(node_options, "N", "Remove nodes covered by fewer than *N* number of path steps.", {'c', "min-depth"});
    args::ValueFlag<uint64_t> max_depth(node_options, "N", "Remove nodes covered by more than *N* number of path steps.", {'C', "max-depth"});
    args::Flag cut_tips(node_options, "bool", "Remove nodes which are graph tips.", {'T', "cut-tips"});
    args::Group edge_opts(parser, "[ Edge Options ]");
    args::Flag edge_depth(edge_opts, "bool", "Remove edges outside of the minimum and maximum coverage rather than"
                                             " nodes. Only set this argument in combination with **-c,"
                                             " –min-coverage**=*N* and **-C, --max-coverage**=*N*.", {'E', "edge-depth"});
    args::ValueFlag<uint64_t> best_edges(edge_opts, "N", "Only keep the *N* most covered inbound and output edges of each node.", {'b', "best-edges"});
    args::Group step_opts(parser, "[ Step Options ]");
    args::ValueFlag<uint64_t> expand_steps(step_opts, "N", "Also include nodes within this many steps of a component passing the prune thresholds.", {'s', "expand-steps"});
    args::ValueFlag<uint64_t> expand_length(step_opts, "N", "Also include nodes within this graph nucleotide distance of a component passing the prune thresholds.", {'l', "expand-length"});
    args::Group path_opts(parser, "[ Path Options ]");
    args::ValueFlag<uint64_t> expand_path_length(path_opts, "N", "Also include nodes within this path length of a component passing the prune thresholds.", {'p', "expand-path-length"});
    args::ValueFlag<std::string> drop_paths(path_opts, "FILE",
                                                  "List of paths to remove. The FILE must "
                                                  "contain one path name per line and a subset of all paths can be specified.",
                                                  {'r', "drop-paths"});
    args::Flag drop_all_paths(path_opts, "bool", "Remove all paths from the graph.", {'D', "drop-all-paths"});
    args::Flag drop_empty_paths(path_opts, "bool", "Remove empty paths from the graph.", {'y', "drop-empty-paths"});
    args::Flag split_paths(path_opts, "bool", "When pruning nodes/edges, keep the surviving paths by splitting them into subpaths (named *PATH:START-END*) at every removed node and edge, instead of dropping all paths.", {'S', "split-paths"});
    args::Flag drop_affected_paths(path_opts, "bool", "When pruning nodes/edges, drop only the paths that traverse a removed node or edge, keeping the untouched paths intact, instead of dropping all paths.", {'A', "drop-affected-paths"});
    args::ValueFlag<uint64_t> cut_tips_min_depth(path_opts, "N", "Remove nodes which are graph tips and have less than *N* path depth.", {'m', "cut-tips-min-depth"});
    args::Flag remove_isolated(path_opts, "bool", "Remove isolated nodes covered by a single path.", {'I', "remove-isolated"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<int> threads(threading_opts, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
	args::Group processing_info_opts(parser, "[ Processing Information ]");
	args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
	args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi prune.", {'h', "help"});

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
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    assert(argc > 0);

    if (!dg_in_file) {
        std::cerr << "[odgi::prune] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (!dg_out_file) {
        std::cerr
                << "[odgi::prune] error: please specify an output file to where to store the pruned graph via -o=[FILE], --out=[FILE]."
                << std::endl;
        return 1;
    }

    // How to treat paths affected by pruning
    const bool split_paths_requested = args::get(split_paths);
    const bool drop_affected_requested = args::get(drop_affected_paths);
    const bool rebuild_paths_requested = split_paths_requested || drop_affected_requested;
    if (rebuild_paths_requested) {
        if (split_paths_requested && drop_affected_requested) {
            std::cerr << "[odgi::prune] error: -S/--split-paths and -A/--drop-affected-paths are mutually exclusive." << std::endl;
            return 1;
        }
        if (!(args::get(max_degree) || args::get(max_furcations) || args::get(min_depth)
              || args::get(max_depth) || args::get(best_edges))) {
            std::cerr << "[odgi::prune] error: -S/--split-paths and -A/--drop-affected-paths require a node/edge pruning filter (-d/--max-degree, -e/--max-furcations, -c/--min-depth, -C/--max-depth, -E/--edge-depth, or -b/--best-edges)." << std::endl;
            return 1;
        }
        if (args::get(expand_path_length)) {
            std::cerr << "[odgi::prune] error: -S/--split-paths and -A/--drop-affected-paths cannot be combined with -p/--expand-path-length, which already preserves paths as segments." << std::endl;
            return 1;
        }
    }

	const int n_threads = threads ? args::get(threads) : 1;

	graph_t graph;
    {
        const std::string infile = args::get(dg_in_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
				utils::handle_gfa_odgi_input(infile, "prune", args::get(progress), n_threads, graph);
            }
        }
    }

    omp_set_num_threads(n_threads);

    // Removing a node/edge breaks the paths through it. By default we drop all paths and warn;
    // with -S/-A we snapshot the paths first (before any filter, since -d clears them), run the
    // filters, then rebuild the survivors edge-aware (algorithms::rebuild_pruned_paths).
    const uint64_t n_paths_before = graph.get_path_count();
    // -c1 alone removes only depth-0 (unpathed) nodes, so it can't damage a path; other filters can.
    const bool only_depth0 = (args::get(min_depth) == 1 && args::get(max_depth) == 0)
        && !args::get(max_degree) && !args::get(max_furcations) && !args::get(best_edges);
    const bool damages_paths = n_paths_before > 0 && !only_depth0
        && (args::get(max_degree) || args::get(max_furcations)
            || args::get(min_depth) || args::get(max_depth) || args::get(best_edges));
    const bool rebuild = rebuild_paths_requested && damages_paths;

    std::vector<algorithms::recorded_path_t> recorded_paths;
    if (rebuild) {
        recorded_paths = algorithms::record_paths(graph);
    }

    if (args::get(max_degree)) {
        graph.clear_paths();
        algorithms::remove_high_degree_nodes(graph, args::get(max_degree));
    }
    if (args::get(max_furcations)) {
        std::vector<edge_t> to_prune = algorithms::find_edges_to_prune(graph, args::get(kmer_length), args::get(max_furcations), n_threads);
        for (auto& edge : to_prune) {
            graph.destroy_edge(edge);
        }
        // leaves edge-damaged paths in place; the path handling below drops or rebuilds them
    }
    if (args::get(min_depth) || args::get(max_depth) || args::get(best_edges)) {
        std::vector<handle_t> handles_to_drop;
        std::vector<edge_t> edges_to_drop_depth;
        std::vector<edge_t> edges_to_drop_best;

        if (args::get(edge_depth)) {
            edges_to_drop_depth = algorithms::find_edges_exceeding_depth_limits(graph, args::get(min_depth), args::get(max_depth));
        } else {
            handles_to_drop = algorithms::find_handles_exceeding_depth_limits(graph, args::get(min_depth), args::get(max_depth));
        }
        if (args::get(best_edges)) {
            edges_to_drop_best = algorithms::keep_mutual_best_edges(graph, args::get(best_edges));
        }
        auto do_destroy =
            [&]() {
                if (!only_depth0) {
                    graph.clear_paths();
                }
                for (auto& edge : edges_to_drop_depth) {
                    graph.destroy_edge(edge);
                }
                for (auto& edge : edges_to_drop_best) {
                    graph.destroy_edge(edge);
                }
                for (auto& handle : handles_to_drop) {
                    graph.destroy_handle(handle);
                }
            };
        // A copy of the graph topology is needed to expand the pruned graph against the original.
        const bool need_source =
            args::get(expand_steps) || args::get(expand_length) || args::get(expand_path_length);
        graph_t source;
        if (need_source) {
            source.copy(graph);
        }
        do_destroy();
        if (args::get(expand_steps)) {
            algorithms::expand_context(&source, &graph, args::get(expand_steps), true);
        } else if (args::get(expand_length)) {
            algorithms::expand_context(&source, &graph, args::get(expand_length), false);
        } else if (args::get(expand_path_length)) {
            // expand_context_with_paths re-embeds the source paths as fresh segments and requires
            // an empty target; do_destroy() leaves the paths in place for the -c1 case, so clear
            // them here to avoid a duplicate path-name error.
            graph.clear_paths();
            algorithms::expand_context_with_paths(&source, &graph, args::get(expand_path_length), false);
        }
    }

    // Rebuild or drop the damaged paths, once, for whichever filters ran above.
    if (rebuild) {
        graph.clear_paths(); // drop any partial/invalid paths a filter left behind (e.g. -e)
        const algorithms::damaged_path_policy_t policy = split_paths_requested
            ? algorithms::damaged_path_policy_t::split
            : algorithms::damaged_path_policy_t::drop_affected;
        const algorithms::prune_path_rebuild_stats_t st =
            algorithms::rebuild_pruned_paths(recorded_paths, graph, policy);
        if (split_paths_requested) {
            std::cerr << "[odgi::prune] split " << st.paths_split << " damaged path(s) into "
                      << st.subpaths_created << " subpath(s), kept " << st.paths_intact
                      << " path(s) intact";
            if (st.paths_dropped) {
                std::cerr << ", dropped " << st.paths_dropped << " fully removed path(s)";
            }
            std::cerr << "." << std::endl;
        } else {
            std::cerr << "[odgi::prune] kept " << st.paths_intact << " intact path(s), dropped "
                      << st.paths_dropped << " affected path(s)." << std::endl;
        }
    } else if (damages_paths && !args::get(expand_path_length)) {
        // Default: ensure the damaged paths are gone (some filters, e.g. -e, keep them) and warn.
        // -p re-embeds paths itself, so it is excluded here.
        if (graph.get_path_count() > 0) {
            graph.clear_paths();
        }
        std::cerr << "[odgi::prune] warning: pruning removed nodes/edges used by the paths, so all "
                  << n_paths_before << " path(s) were dropped from the output. "
                  << "Re-run with -S/--split-paths (keep them as subpaths) or "
                  << "-A/--drop-affected-paths (drop only the touched paths)." << std::endl;
    }
    if (args::get(cut_tips)) {
        algorithms::cut_tips(graph, args::get(cut_tips_min_depth));
        graph.optimize();
    }
    if (args::get(remove_isolated)) {
        algorithms::remove_isolated_paths(graph);
        graph.optimize();
    }
    if (args::get(drop_all_paths)) {
        graph.clear_paths();
    }
    if (args::get(drop_empty_paths)) {
        std::vector<path_handle_t> empty_paths;
        empty_paths.reserve(graph.get_path_count());
        graph.for_each_path_handle([&](const path_handle_t& path) {
            if (graph.is_empty(path)) {
                empty_paths.push_back(path);
            }
        });

        for (auto path: empty_paths) {
            graph.destroy_path(path);
        }
    } else if (!args::get(drop_paths).empty()){
        std::vector<path_handle_t> paths_to_remove;

        std::ifstream path_names_in(args::get(drop_paths));

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
                        paths_to_remove.push_back(path);
                    } else {
                        std::cerr << "[odgi::prune] error: in the path list there are duplicated path names."
                                  << std::endl;
                        exit(1);
                    }
                }

                ++num_of_paths_in_file;
            }
        }

        path_names_in.close();

        std::cerr << "[odgi::prune] found " << paths_to_remove.size() << "/" << num_of_paths_in_file
                  << " paths to remove." << std::endl;

        for(auto& path : paths_to_remove) {
            graph.destroy_path(path);
        }
    }

    {
        const std::string outfile = args::get(dg_out_file);
        if (!outfile.empty()) {
            if (outfile == "-") {
                graph.serialize(std::cout);
            } else {
                ofstream f(outfile.c_str());
                graph.serialize(f);
                f.close();
            }
        }
    }

    return 0;
}

static Subcommand odgi_prune("prune", "Remove parts of the graph.",
                              PIPELINE, 3, main_prune);


}
