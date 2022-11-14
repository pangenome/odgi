#include <xp.hpp>
#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/topological_sort.hpp"
#include "algorithms/cycle_breaking_sort.hpp"
#include "algorithms/id_ordered_paths.hpp"
#include "algorithms/dagify.hpp"
#include "algorithms/split_strands.hpp"
#include "algorithms/dagify_sort.hpp"
#include "algorithms/random_order.hpp"
//#include "algorithms/mondriaan_sort.hpp"
#include "algorithms/xp.hpp"
#include "algorithms/pg_sgd/path_sgd.hpp"
#include "algorithms/pg_sgd/path_sgd_helper.hpp"
#include "algorithms/pg_sgd/path_sgd_ssi.hpp"
#include "algorithms/groom.hpp"
#include "algorithms/stepindex.hpp"
#include "utils.hpp"

namespace odgi {

using namespace odgi::subcommand;

#define MAX_NUMBER_OF_ZIPF_DISTRIBUTIONS 100

int main_sort(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi sort";
    argv[0] = (char*)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("Apply different kind of sorting algorithms to a graph. The most prominent one is the PG-SGD sorting algorithm.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::ValueFlag<std::string> dg_out_file(mandatory_opts, "FILE", "Write the sorted dynamic succinct variation graph to this file. A file"
                                                             " ending with *.og* is recommended.", {'o', "out"});
    args::Group files_io_opts(parser, "[ Files IO Options ]");
    args::ValueFlag<std::string> xp_in_file(files_io_opts, "FILE", "Load the succinct variation graph index from this *FILE*. The file name usually ends with *.xp*.", {'X', "path-index"});
	/// @Andrea new argument
	args::ValueFlag<std::string> ssi_in_file(files_io_opts, "FILE", "Load the sampled step index from this *FILE*. The file name usually ends with *.ssi*.", {'e', "sampled-step-index"});
    args::ValueFlag<std::string> sort_order_in(files_io_opts, "FILE", "*FILE* containing the sort order. Each line contains one node identifer.", {'s', "sort-order"});
    args::ValueFlag<std::string> tmp_base(files_io_opts, "PATH", "directory for temporary files", {'C', "temp-dir"});
    args::Group topo_sorts_opts(parser, "[ Topological Sort Options ]");
    args::Flag breadth_first(topo_sorts_opts, "breadth_first", "Use a (chunked) breadth first topological sort.", {'b', "breadth-first"});
    args::ValueFlag<uint64_t> breadth_first_chunk(topo_sorts_opts, "N", "Chunk size for breadth first topological sort. Specify how many"
                                                                        " nucleotides to grap at once in each BFS phase.", {'B', "breadth-first-chunk"});
    args::Flag cycle_breaking(topo_sorts_opts, "cycle_breaking", "Use a cycle breaking sort.", {'c', "cycle-breaking"});
    args::Flag depth_first(topo_sorts_opts, "depth_first", "Use a (chunked) depth first topological sort.", {'z', "depth-first"});
    args::ValueFlag<uint64_t> depth_first_chunk(topo_sorts_opts, "N", "Chunk size for the depth first topological sort. Specify how many"
                                                                      " nucleotides to grap at once in each DFS phase.", {'Z', "depth-first-chunk"});
    args::Flag two(topo_sorts_opts, "two", "Use a two-way topological algorithm for sorting. It is a maximum of"
                                           " head-first and tail-first topological sort.", {'w', "two-way"});
    args::Flag no_seeds(topo_sorts_opts, "no-seeds", "Don't use heads or tails to seed the topological sort.", {'n', "no-seeds"});
    // other sorts
    args::Group random_sort_opts(parser, "[ Random Sort Options ]");
    args::Flag randomize(random_sort_opts, "random", "Randomly sort the graph.", {'r', "random"});
    args::Group dagify_sort_opts(parser, "[ DAGify Sort Options ]");
    args::Flag dagify(dagify_sort_opts, "dagify", "Sort on the basis of a DAGified graph.", {'d', "dagify-sort"});
    /// path guided linear 1D SGD
    args::Group pg_sgd_opts(parser, "[ Path Guided 1D SGD Sort ]");
    args::Flag p_sgd(pg_sgd_opts, "path-sgd", "Apply the path-guided linear 1D SGD algorithm to organize graph.", {'Y', "path-sgd"});
	/// @Andrea new argument
	args::Flag p_sgd_xp(pg_sgd_opts, "path-sgd", "Use the faster but memory expensive XP index.", {'m', "path-sgd-path-index"});
    args::ValueFlag<std::string> p_sgd_in_file(pg_sgd_opts, "FILE", "Specify a line separated list of paths to sample from for the on the"
                                                                    " fly term generation process in the path guided linear 1D SGD (default: sample from all paths).", {'f', "path-sgd-use-paths"});
    args::ValueFlag<double> p_sgd_min_term_updates_paths(pg_sgd_opts, "N", "The minimum number of terms to be updated before a new path guided"
                                                                           " linear 1D SGD iteration with adjusted learning rate eta starts,"
                                                                           " expressed as a multiple of total path steps (default: *1.0*). Can be overwritten by *-U, -path-sgd-min-term-updates-nodes=N*.", {'G', "path-sgd-min-term-updates-paths"});
    args::ValueFlag<double> p_sgd_min_term_updates_num_nodes(pg_sgd_opts, "N", "The minimum number of terms to be updated before a new path guided"
                                                                               " linear 1D SGD iteration with adjusted learning rate eta starts,"
                                                                               " expressed as a multiple of the number of nodes (default: NONE. *-G,path-sgd-min-term-updates-paths=N* is used).", {'U', "path-sgd-min-term-updates-nodes"});
    args::ValueFlag<double> p_sgd_delta(pg_sgd_opts, "N", "The threshold of maximum displacement approximately in bp at which to"
                                                          " stop path guided linear 1D SGD (default: *0.0*).", {'j', "path-sgd-delta"});
    args::ValueFlag<double> p_sgd_eps(pg_sgd_opts, "N", "The final learning rate for path guided linear 1D SGD model (default: *0.01*).", {'g', "path-sgd-eps"});
    args::ValueFlag<double> p_sgd_eta_max(pg_sgd_opts, "N", "The first and maximum learning rate for path guided linear 1D SGD"
                                                            " model (default: *squared steps of longest path in graph*).", {'v', "path-sgd-eta-max"});
    args::ValueFlag<double> p_sgd_zipf_theta(pg_sgd_opts, "N", "The theta value for the Zipfian distribution which is used as the"
                                                               " sampling method for the second node of one term in the path guided"
                                                               " linear 1D SGD model (default: *0.99*).", {'a', "path-sgd-zipf-theta"});
    args::ValueFlag<uint64_t> p_sgd_iter_max(pg_sgd_opts, "N", "The maximum number of iterations for path guided linear 1D SGD model (default: 100).", {'x', "path-sgd-iter-max"});
    args::ValueFlag<double> p_sgd_cooling(pg_sgd_opts, "N",
                                          "Use this fraction of the iterations for layout annealing (default: 0.5).",
                                          {'K', "path-sgd-cooling"});
    args::ValueFlag<uint64_t> p_sgd_iter_with_max_learning_rate(pg_sgd_opts, "N", "The iteration where the learning rate is max for path guided linear 1D SGD model (default: *0*).", {'F', "iteration-max-learning-rate"});
    args::ValueFlag<uint64_t> p_sgd_zipf_space(pg_sgd_opts, "N", "The maximum space size of the Zipfian distribution which is used as"
                                                                 " the sampling method for the second node of one term in the path guided"
                                                                 " linear 1D SGD model (default: *longest path length*).", {'k', "path-sgd-zipf-space"});
    args::ValueFlag<uint64_t> p_sgd_zipf_space_max(pg_sgd_opts, "N", "The maximum space size of the Zipfian distribution beyond which"
                                                                     " quantization occurs (default: *100*).", {'I', "path-sgd-zipf-space-max"});
    args::ValueFlag<uint64_t> p_sgd_zipf_space_quantization_step(pg_sgd_opts, "N", "Quantization step size when the maximum space size of the Zipfian"
                                                                                   " distribution is exceeded (default: *100*).", {'l', "path-sgd-zipf-space-quantization-step"});
    args::ValueFlag<uint64_t> p_sgd_zipf_max_number_of_distributions(pg_sgd_opts, "N", "Approximate maximum number of Zipfian distributions to calculate (default: *100*).", {'y', "path-sgd-zipf-max-num-distributions"});
    args::ValueFlag<std::string> p_sgd_seed(pg_sgd_opts, "STRING", "| Set the seed for the deterministic 1-threaded path guided linear 1D SGD model (default: *pangenomic!*).", {'q', "path-sgd-seed"});
    args::ValueFlag<std::string> p_sgd_snapshot(pg_sgd_opts, "STRING", "Set the prefix to which each snapshot graph of a path guided 1D SGD"
                                                                       " iteration should be written to. This is turned off per default. This"
                                                                       " argument only works when *-Y, –path-sgd* was specified. Not applicable"
                                                                       " in a pipeline of sorts.", {'u', "path-sgd-snapshot"});
	args::ValueFlag<std::string> _p_sgd_target_paths(pg_sgd_opts, "FILE", "Read the paths that should be considered as target paths (references) from this *FILE*. PG-SGD will keep the nodes of the given paths fixed. A path's rank determines it's weight for decision making and is given by its position in the given *FILE*.", {'H', "target-paths"});

	/// pipeline
    args::Group pipeline_sort_opts(parser, "[ Pipeline Sorting Options ]");
    args::ValueFlag<std::string> pipeline(pipeline_sort_opts, "STRING", "Apply a series of sorts, based on single character command line"
                                                                        " arguments given to this command (default: NONE). *s*: Topolocigal sort, heads only. *n*: Topological sort, no heads, no tails. *d*: DAGify sort. *c*: Cycle breaking sort. *b*: Breadth first topological sort. *z*: Depth first topological sort. *w*: Two-way topological sort. *r*: Random sort. *Y*: PG-SGD 1D sort. *f*: Reverse order. *g*: Groom the graph. An example could be *Ygs*.", {'p', "pipeline"});
    /// paths
    args::Group path_sorting_opts(parser, "[ Path Sorting Options ]");
    args::Flag paths_by_min_node_id(path_sorting_opts, "paths-min", "Sort paths by their lowest contained node identifier.", {'L', "paths-min"});
    args::Flag paths_by_max_node_id(path_sorting_opts, "paths-max", "Sort paths by their highest contained node identifier.", {'M', "paths-max"});
    args::Flag paths_by_avg_node_id(path_sorting_opts, "paths-avg", "Sort paths by their average contained node identifier.", {'A', "paths-avg"});
    args::Flag paths_by_avg_node_id_rev(path_sorting_opts, "paths-avg-rev", "Sort paths in reverse by their average contained node identifier.", {'R', "paths-avg-rev"});
    args::ValueFlag<std::string> path_delim(path_sorting_opts, "path-delim", "Sort paths in bins by their prefix up to this delimiter.", {'D', "path-delim"});
    /// misc
    args::Group optimize_opts(parser, "[ Optimize Options ]");
    args::Flag optimize(optimize_opts, "optimize", "Use the MutableHandleGraph::optimize method to compact the node"
                                                   " identifier space.", {'O', "optimize"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> nthreads(threading_opts, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
    args::Group processing_info_opts(parser, "[ Processing Information ]");
    args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi sort.", {'h', "help"});

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

    if (!dg_in_file) {
        std::cerr << "[odgi::sort] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (!dg_out_file || args::get(dg_out_file).empty()) {
        std::cerr << "[odgi::sort] error: please specify an output file to store the graph via -o=[FILE], --out=[FILE]." << std::endl;
        return 1;
    }

    if (args::get(p_sgd_zipf_space_quantization_step) && args::get(p_sgd_zipf_max_number_of_distributions)){
        std::cerr
                << "[odgi::sort] error: please specify -l/--path-sgd-zipf-space-quantization-step or -y/--path-sgd-zipf-max-num-distributions, not both."
                << std::endl;
        return 1;
    }

	const bool xp_way = ((xp_in_file || p_sgd_xp) && p_sgd) ;

	if (xp_way && ssi_in_file) {
		std::cerr
				<< "[odgi::sort] error: please specify exclusively -X=[FILE]/--path-index=[FILE], -m/--path-sgd-path-index, or -e=[FILE]/--sampled-step-index=[FILE]."
				<< std::endl;
		return 1;
	}

	if (xp_in_file && p_sgd_xp) {
		std::cerr
				<< "[odgi::sort] error: please specify -X=[FILE]/--path-index=[FILE] or -m=/--path-sgd-path-index, not both."
				<< std::endl;
		return 1;
	}

	const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

	graph_t graph;
    assert(argc > 0);
    std::string infile = args::get(dg_in_file);
    if (infile.size()) {
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
			utils::handle_gfa_odgi_input(infile, "sort", args::get(progress), num_threads, graph);
        }
    }

    // default settings
    uint64_t df_chunk_size = args::get(depth_first_chunk) ? args::get(depth_first_chunk) : 1000;
    uint64_t bf_chunk_size = args::get(breadth_first_chunk) ? args::get(breadth_first_chunk) : std::numeric_limits<uint64_t>::max();

    graph.set_number_of_threads(num_threads);

    // default parameters
    std::string path_sgd_seed;
    if (p_sgd_seed) {
        if (num_threads > 1) {
            std::cerr << "[odgi::sort] error: please only specify a seed for the path guided 1D linear SGD when using 1 thread." << std::endl;
            return 1;
        }
        path_sgd_seed = args::get(p_sgd_seed);
    } else {
        path_sgd_seed = "pangenomic!";
    }
    if (p_sgd_min_term_updates_paths && p_sgd_min_term_updates_num_nodes) {
        std::cerr << "[odgi::sort] error: there can only be one argument provided for the minimum number of term updates in the path guided 1D SGD."
                     "Please either use -G=[N], path-sgd-min-term-updates-paths=[N] or -U=[N], path-sgd-min-term-updates-nodes=[N]." << std::endl;
        return 1;
    }

	// SSI should be the new default, else we take the XP we received from the input or if the input is an empty string we generate a new one anyhow
	// if an XP is given or we want to work with an XP based PG-SGD, we should set the temporary file
	if (xp_way) {
		std::cerr << "[odgi::sort] using XP index for the PG-SGD algorithm" << std::endl;
		if (tmp_base) {
			xp::temp_file::set_dir(args::get(tmp_base));
		} else {
			char* cwd = get_current_dir_name();
			xp::temp_file::set_dir(std::string(cwd));
			free(cwd);
		}
	}

    // If required, first of all, optimize the graph so that it is optimized for subsequent algorithms (if required)
    if (args::get(optimize)) {
        graph.optimize();
    }

    uint64_t path_sgd_iter_max = args::get(p_sgd_iter_max) ? args::get(p_sgd_iter_max) : 100;
    uint64_t path_sgd_iter_max_learning_rate = args::get(p_sgd_iter_with_max_learning_rate) ? args::get(p_sgd_iter_with_max_learning_rate) : 0;
    double path_sgd_zipf_theta = args::get(p_sgd_zipf_theta) ? args::get(p_sgd_zipf_theta) : 0.99;
    double path_sgd_eps = args::get(p_sgd_eps) ? args::get(p_sgd_eps) : 0.01;
    double path_sgd_delta = args::get(p_sgd_delta) ? args::get(p_sgd_delta) : 0;
    double path_sgd_max_eta = 0; // update below
    //double path_sgd_cooling_start = 2.0; // disabled
    double path_sgd_cooling = p_sgd_cooling ? args::get(p_sgd_cooling) : 0.5;
    // will be filled, if the user decides to write a snapshot of the graph after each sorting iteration
    std::vector<std::string> snapshots;
    const bool snapshot = p_sgd_snapshot;
    // default parameters that need a path index to be present
    uint64_t path_sgd_min_term_updates;
    uint64_t path_sgd_zipf_space, path_sgd_zipf_space_max, path_sgd_zipf_space_quantization_step, path_sgd_zipf_max_number_of_distributions;
    std::vector<path_handle_t> path_sgd_use_paths;
    xp::XP path_index;
	algorithms::step_index_t sampled_step_index;
    bool fresh_step_index = false;
    std::string snapshot_prefix;
    if (snapshot) {
        snapshot_prefix = args::get(p_sgd_snapshot);
    }
	// do we only want so sample from a subset of paths?
	if (p_sgd_in_file) {
		std::string buf;
		std::ifstream use_paths(args::get(p_sgd_in_file).c_str());
		while (std::getline(use_paths, buf)) {
			// check if the path is actually in the graph, else print an error and exit 1
			if (graph.has_path(buf)) {
				path_sgd_use_paths.push_back(graph.get_path_handle(buf));
			} else {
				std::cerr << "[odgi::sort] error: path '" << buf
						  << "' as was given by -f=[FILE], --path-sgd-use-paths=[FILE]"
							 " is not present in the graph. Please remove this path from the file and restart 'odgi sort'.";
			}
		}
		use_paths.close();
	} else {
		graph.for_each_path_handle(
				[&](const path_handle_t &path) {
					path_sgd_use_paths.push_back(path);
				});
	}
	std::vector<bool> is_ref;
	std::vector<path_handle_t> target_paths;
	/// @Andrea
	uint64_t ssi_sample_rate = 8;
	uint64_t sum_path_step_count;
    if (p_sgd || args::get(pipeline).find('Y') != std::string::npos) {
		if (_p_sgd_target_paths) {
			target_paths = utils::load_paths(args::get(_p_sgd_target_paths), graph, "sort");
			algorithms::sort_graph_by_target_paths(graph, target_paths, is_ref, args::get(progress));
		}
		// SSI should be the new default, else we take the XP we received from the input or if the input is an empty string we generate a new one anyhow
		/// @Andrea
        if (xp_in_file) {
            std::ifstream in;
            in.open(args::get(xp_in_file));
            path_index.load(in);
            in.close();
        } else if (p_sgd_xp) {
			// only if we set the flag!
            path_index.from_handle_graph(graph, num_threads);
			// do we load the SSI from file?
        } else if (ssi_in_file) {
			// FIXME we could also only allow an integer as input to build an ssi for the given sample rate? @Andrea
			sampled_step_index.load(args::get(ssi_in_file));
			// store the SSI sampling rate
			ssi_sample_rate = sampled_step_index.get_sample_rate();
			// or do we build it from scratch?
		} else {
			// stepindex only for selected paths
			sampled_step_index.from_handle_graph(graph, path_sgd_use_paths, num_threads, args::get(progress), ssi_sample_rate);
		}
        fresh_step_index = true;
		/// @Andrea
        sum_path_step_count = algorithms::get_sum_path_step_count(path_sgd_use_paths, graph);
        if (args::get(p_sgd_min_term_updates_paths)) {
            path_sgd_min_term_updates = args::get(p_sgd_min_term_updates_paths) * sum_path_step_count;
        } else {
            if (args::get(p_sgd_min_term_updates_num_nodes)) {
                path_sgd_min_term_updates = args::get(p_sgd_min_term_updates_num_nodes) * graph.get_node_count();
            } else {
                path_sgd_min_term_updates = 1.0 * sum_path_step_count;
            }
        }
		/// @Andrea
        uint64_t max_path_step_count = algorithms::get_max_path_step_count(path_sgd_use_paths, graph);
		// we might have to use the SSI here
		if (xp_way) {
			path_sgd_zipf_space = args::get(p_sgd_zipf_space) ? args::get(p_sgd_zipf_space) : algorithms::get_max_path_length_xp(path_sgd_use_paths, path_index);
		} else {
			path_sgd_zipf_space = args::get(p_sgd_zipf_space) ? args::get(p_sgd_zipf_space) : algorithms::get_max_path_length_ssi(path_sgd_use_paths, sampled_step_index);
		}
        path_sgd_zipf_space_max = args::get(p_sgd_zipf_space_max) ? args::get(p_sgd_zipf_space_max) : 100;

        path_sgd_zipf_max_number_of_distributions = args::get(p_sgd_zipf_max_number_of_distributions) ? std::max(
                (uint64_t) path_sgd_zipf_space_max + 1,
                (uint64_t) args::get(p_sgd_zipf_max_number_of_distributions)
        ) : std::max((uint64_t) path_sgd_zipf_space_max + 1,
                     (uint64_t) MAX_NUMBER_OF_ZIPF_DISTRIBUTIONS);

        if (args::get(progress)) {
            std::cerr << "path_sgd_zipf_space_max: " << path_sgd_zipf_space_max << std::endl;
            std::cerr << "path_sgd_zipf_max_number_of_distributions: " << path_sgd_zipf_max_number_of_distributions << std::endl;
        }

        if (args::get(p_sgd_zipf_space_quantization_step)) {
            path_sgd_zipf_space_quantization_step = std::max((uint64_t) 2, args::get(p_sgd_zipf_space_quantization_step));
        } else {
            if (path_sgd_zipf_space > path_sgd_zipf_space_max && path_sgd_zipf_max_number_of_distributions > path_sgd_zipf_space_max) {
                path_sgd_zipf_space_quantization_step = std::max(
                        (uint64_t) 2,
                        (uint64_t) ceil( (double) (path_sgd_zipf_space - path_sgd_zipf_space_max) / (double) (path_sgd_zipf_max_number_of_distributions - path_sgd_zipf_space_max))
                );
            } else {
                path_sgd_zipf_space_quantization_step = 100;
            }
        }

        path_sgd_max_eta = args::get(p_sgd_eta_max) ? args::get(p_sgd_eta_max) : max_path_step_count * max_path_step_count;
    }

    // did we groom the graph?
    if (!args::get(pipeline).empty()) {
        // for each sort type, apply it to the graph
        std::vector<handle_t> order;
        for (auto c : args::get(pipeline)) {
            switch (c) {
                case 's':
                    order = algorithms::topological_order(&graph, true, false, args::get(progress));
                    break;
                case 'n':
                    order = algorithms::topological_order(&graph, false, false, args::get(progress));
                    break;
                case 'd': {
                    graph_t split, into;
                    order = algorithms::dagify_sort(graph, split, into);
                }
                    break;
                case 'c':
                    order = algorithms::cycle_breaking_sort(graph);
                    break;
                case 'b':
                    order = algorithms::breadth_first_topological_order(graph, bf_chunk_size);
                    break;
                case 'z':
                    order = algorithms::depth_first_topological_order(graph, df_chunk_size);
                    break;
                case 'w':
                    order = algorithms::two_way_topological_order(&graph);
                    break;
                case 'r':
                    order = algorithms::random_order(graph);
                    break;
                case 'Y': {
                    if (!fresh_step_index) {
						if (xp_way) {
							path_index.clean();
						} else {
							sampled_step_index.clean();
						}
						// do we have to sort by reference nodes first?
						if (_p_sgd_target_paths) {
							algorithms::sort_graph_by_target_paths(graph, target_paths, is_ref, args::get(progress));
						}
						if (xp_way) {
							path_index.from_handle_graph(graph, num_threads);
						} else {
							sampled_step_index.from_handle_graph(graph, path_sgd_use_paths, num_threads, args::get(progress), ssi_sample_rate);
						}
                    }
					if (xp_way) {
						order = algorithms::path_linear_sgd_order(graph,
																  path_index,
																  path_sgd_use_paths,
																  path_sgd_iter_max,
																  path_sgd_iter_max_learning_rate,
																  path_sgd_min_term_updates,
																  path_sgd_delta,
																  path_sgd_eps,
																  path_sgd_max_eta,
																  path_sgd_zipf_theta,
																  path_sgd_zipf_space,
																  path_sgd_zipf_space_max,
																  path_sgd_zipf_space_quantization_step,
																  path_sgd_cooling,
																  num_threads,
																  progress,
																  path_sgd_seed,
																  snapshot,
																  snapshot_prefix,
																  _p_sgd_target_paths,
																  is_ref);
					} else {
						order = algorithms::path_linear_sgd_order_ssi(graph,
																  sampled_step_index,
																  path_sgd_use_paths,
																  path_sgd_iter_max,
																  path_sgd_iter_max_learning_rate,
																  path_sgd_min_term_updates,
																  path_sgd_delta,
																  path_sgd_eps,
																  path_sgd_max_eta,
																  path_sgd_zipf_theta,
																  path_sgd_zipf_space,
																  path_sgd_zipf_space_max,
																  path_sgd_zipf_space_quantization_step,
																  path_sgd_cooling,
																  num_threads,
																  progress,
																  path_sgd_seed,
																  snapshot,
																  snapshot_prefix,
																  _p_sgd_target_paths,
																  is_ref,
																  sum_path_step_count);
					}
                    break;
                }
                case 'f':
                    order.clear();
                    graph.for_each_handle([&order](const handle_t &handle) {
                        order.push_back(handle);
                    });
                    std::reverse(order.begin(), order.end());
                    break;
                case 'g': {
                    order = algorithms::groom(graph, progress, target_paths);
                    break;
                }
                default:
                    break;
            }
            if (order.size() != graph.get_node_count()) {
                std::cerr << "[odgi::sort] error: expected " << graph.get_node_count()
                          << " handles in the order "
                          << "but got " << order.size() << std::endl;
                assert(false);
            }
            graph.apply_ordering(order, true);
            fresh_step_index = false;
        }
    } else if (args::get(two)) {
        graph.apply_ordering(algorithms::two_way_topological_order(&graph), true);
    } else if (!args::get(sort_order_in).empty()) {
        std::vector<handle_t> given_order;
        std::string buf;
        std::ifstream in_order(args::get(sort_order_in).c_str());
        while (std::getline(in_order, buf)) {
            given_order.push_back(graph.get_handle(std::stol(buf)));
        }
        graph.apply_ordering(given_order, true);
    } else if (args::get(dagify)) {
        graph_t split, into;
        graph.apply_ordering(algorithms::dagify_sort(graph, split, into), true);
    } else if (args::get(cycle_breaking)) {
        graph.apply_ordering(algorithms::cycle_breaking_sort(graph), true);
    } else if (args::get(no_seeds)) {
        graph.apply_ordering(algorithms::topological_order(&graph, false, false, args::get(progress)), true);
    } else if (args::get(p_sgd)) {
		// FIXME SSI should be the new default, else we take the XP we received from the input or if the input is an empty string we generate a new one anyhow
		std::vector<handle_t> order;
		if (xp_way) {
			order = algorithms::path_linear_sgd_order(graph,
													  path_index,
													  path_sgd_use_paths,
													  path_sgd_iter_max,
													  path_sgd_iter_max_learning_rate,
													  path_sgd_min_term_updates,
													  path_sgd_delta,
													  path_sgd_eps,
													  path_sgd_max_eta,
													  path_sgd_zipf_theta,
													  path_sgd_zipf_space,
													  path_sgd_zipf_space_max,
													  path_sgd_zipf_space_quantization_step,
													  path_sgd_cooling,
													  num_threads,
													  progress,
													  path_sgd_seed,
													  snapshot,
													  snapshot_prefix,
													  _p_sgd_target_paths,
													  is_ref);
		} else {
			order = algorithms::path_linear_sgd_order_ssi(graph,
													  sampled_step_index,
													  path_sgd_use_paths,
													  path_sgd_iter_max,
													  path_sgd_iter_max_learning_rate,
													  path_sgd_min_term_updates,
													  path_sgd_delta,
													  path_sgd_eps,
													  path_sgd_max_eta,
													  path_sgd_zipf_theta,
													  path_sgd_zipf_space,
													  path_sgd_zipf_space_max,
													  path_sgd_zipf_space_quantization_step,
													  path_sgd_cooling,
													  num_threads,
													  progress,
													  path_sgd_seed,
													  snapshot,
													  snapshot_prefix,
													  _p_sgd_target_paths,
													  is_ref,
													  sum_path_step_count);
		}
        graph.apply_ordering(order, true);
    } else if (args::get(breadth_first)) {
        graph.apply_ordering(algorithms::breadth_first_topological_order(graph, bf_chunk_size), true);
    } else if (args::get(depth_first)) {
        graph.apply_ordering(algorithms::depth_first_topological_order(graph, df_chunk_size), true);
    } else if (args::get(randomize)) {
        graph.apply_ordering(algorithms::random_order(graph), true);
    } else {
        // To be able to only optimize the graph, avoiding the topological sorting if nothing else is requested
        if (!args::get(optimize)) {
            graph.apply_ordering(algorithms::topological_order(&graph, true, false, args::get(progress)), true);
        }
    }
    if (args::get(paths_by_min_node_id)) {
        graph.apply_path_ordering(
                algorithms::prefix_and_id_ordered_paths(graph, args::get(path_delim), false, false));
    }
    if (args::get(paths_by_max_node_id)) {
        graph.apply_path_ordering(
                algorithms::prefix_and_id_ordered_paths(graph, args::get(path_delim), false, true));
    }
    if (args::get(paths_by_avg_node_id)) {
        graph.apply_path_ordering(
                algorithms::prefix_and_id_ordered_paths(graph, args::get(path_delim), true, false));
    }
    if (args::get(paths_by_avg_node_id_rev)) {
        graph.apply_path_ordering(
                algorithms::prefix_and_id_ordered_paths(graph, args::get(path_delim), true, true));
    }
    const std::string outfile = args::get(dg_out_file);
    if (outfile == "-") {
        graph.serialize(std::cout);
    } else {
        ofstream f(outfile.c_str());
        graph.serialize(f);
        f.close();
    }
    return 0;
}

static Subcommand odgi_sort("sort", "Apply different kind of sorting algorithms to a graph.",
                              PIPELINE, 3, main_sort);


}
