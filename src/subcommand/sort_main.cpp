#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/topological_sort.hpp"
#include "algorithms/eades_algorithm.hpp"
#include "algorithms/cycle_breaking_sort.hpp"
#include "algorithms/id_ordered_paths.hpp"
#include "algorithms/dagify.hpp"
#include "algorithms/split_strands.hpp"
#include "algorithms/dagify_sort.hpp"
#include "algorithms/random_order.hpp"
#include "algorithms/mondriaan_sort.hpp"
#include "algorithms/linear_sgd.hpp"
#include "algorithms/xp.hpp"
#include "algorithms/path_sgd.hpp"
#include "algorithms/groom.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_sort(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi sort";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("sort a variation graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the graph in this file", {'o', "out"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> xp_in_file(parser, "FILE", "load the path index from this file", {'X', "path-index"});
    //args::Flag show_sort(parser, "show", "write the sort order mapping", {'S', "show"});
    args::ValueFlag<std::string> sort_order_in(parser, "FILE", "load the sort order from this file", {'s', "sort-order"});
    args::Flag cycle_breaking(parser, "cycle_breaking", "use a cycle breaking sort", {'c', "cycle-breaking"});
    args::Flag breadth_first(parser, "breadth_first", "use a breadth first topological sort", {'b', "breadth-first"});
    args::Flag depth_first(parser, "depth_first", "use a chunked depth first topological sort", {'z', "depth-first"});
    args::ValueFlag<uint64_t> breadth_first_chunk(parser, "N", "chunk size for breadth first topological sort", {'B', "breadth-first-chunk"});
    args::ValueFlag<uint64_t> depth_first_chunk(parser, "N", "chunk size for depth first topological sort", {'Z', "depth-first-chunk"});
    args::Flag dagify(parser, "dagify", "sort on the basis of the DAGified graph", {'d', "dagify-sort"});
    args::Flag eades(parser, "eades", "use eades algorithm", {'e', "eades"});
    /// linear SGD
    args::Flag lsgd(parser, "linear-sgd", "apply 1D (linear) SGD algorithm to organize graph", {'S', "linear-sgd"});
    args::ValueFlag<uint64_t> lsgd_bandwidth(parser, "sgd-bandwidth", "bandwidth of linear SGD model (default: 1000)", {'H', "sgd-bandwidth"});
    args::ValueFlag<double> lsgd_sampling_rate(parser, "sgd-sampling-rate", "sample pairs of nodes with probability distance between them divided by the sampling rate (default: 20)", {'Q', "sgd-sampling-rate"});
    args::Flag lsgd_use_paths(parser, "sgd-use-paths", "use paths to structure internode distances in SGD", {'K', "sgd-use-paths"});
    args::ValueFlag<uint64_t> lsgd_iter_max(parser, "sgd-iter-max", "max number of iterations for linear SGD model (default: 30)", {'T', "sgd-iter-max"});
    args::ValueFlag<double> lsgd_eps(parser, "sgd-eps", "final learning rate for linear SGD model (default: 0.01)", {'V', "sgd-eps"});
    args::ValueFlag<double> lsgd_delta(parser, "sgd-delta", "threshold of maximum node displacement (approximately in bp) at which to stop SGD (default: 0)", {'C', "sgd-delta"});
    args::Flag two(parser, "two", "use two-way (max of head-first and tail-first) topological algorithm", {'w', "two-way"});
    args::Flag randomize(parser, "random", "randomly sort the graph", {'r', "random"});
    args::Flag no_seeds(parser, "no-seeds", "don't use heads or tails to seed topological sort", {'n', "no-seeds"});
    args::Flag mondriaan(parser, "mondriaan", "use sparse matrix diagonalization to sort the graph", {'m', "mondriaan"});
    args::ValueFlag<uint64_t> mondriaan_n_parts(parser, "N", "number of partitions for mondriaan", {'N', "mondriaan-n-parts"});
    args::ValueFlag<double> mondriaan_epsilon(parser, "N", "epsilon parameter to mondriaan", {'E', "mondriaan-epsilon"});
    args::Flag mondriaan_path_weight(parser, "path-weight", "weight mondriaan input matrix by path coverage of edges", {'W', "mondriaan-path-weight"});
    /// path guided linear 1D SGD
    args::Flag p_sgd(parser, "path-sgd", "apply path guided linear 1D SGD algorithm to organize graph", {'Y', "path-sgd"});
    args::Flag p_sgd_sample_from_nodes(parser, "path-sgd-sample-from-nodes", "instead of sampling the first term from all nucleotide positions of the paths we sample from all nodes of the paths (default: flag not set)", {'J', "path-sgd-sample-from-nodes"});
    args::ValueFlag<std::string> p_sgd_in_file(parser, "FILE", "specify a line separated list of paths to sample from for the on the fly term generation process in the path guided linear 1D SGD (default: sample from all paths)", {'f', "path-sgd-use-paths"});
    args::ValueFlag<double> p_sgd_min_term_updates_paths(parser, "N", "minimum number of terms to be updated before a new path guided linear 1D SGD iteration with adjusted learning rate eta starts, expressed as a multiple of total path length (default: 0.1)", {'G', "path-sgd-min-term-updates-paths"});
    args::ValueFlag<double> p_sgd_min_term_updates_num_nodes(parser, "N", "minimum number of terms to be updated before a new path guided linear 1D SGD iteration with adjusted learning rate eta starts, expressed as a multiple of the number of nodes (default: argument is not set, the default of -G=[N], path-sgd-min-term-updates-paths=[N] is used)", {'U', "path-sgd-min-term-updates-nodes"});
    args::ValueFlag<double> p_sgd_delta(parser, "N", "threshold of maximum displacement approximately in bp at which to stop path guided linear 1D SGD (default: 0)", {'j', "path-sgd-delta"});
    args::ValueFlag<double> p_sgd_eps(parser, "N", "final learning rate for path guided linear 1D SGD model (default: 0.01)", {'g', "path-sgd-eps"});
    args::ValueFlag<double> p_sgd_eta_max(parser, "N", "first and maximum learning rate for path guided linear 1D SGD model (default: number of nodes in the graph)", {'v', "path-sgd-eta-max"});
    args::ValueFlag<double> p_sgd_zipf_theta(parser, "N", "the theta value for the Zipfian distribution which is used as the sampling method for the second node of one term in the path guided linear 1D SGD model (default: 0.99)", {'a', "path-sgd-zipf-theta"});
    args::ValueFlag<uint64_t> p_sgd_iter_max(parser, "N", "max number of iterations for path guided linear 1D SGD model (default: 30)", {'x', "path-sgd-iter-max"});
    args::ValueFlag<uint64_t> p_sgd_iter_with_max_learning_rate(parser, "N", "iteration where the learning rate is max for path guided linear 1D SGD model (default: 0)", {'F', "iteration-max-learning-rate"});
    args::ValueFlag<uint64_t> p_sgd_zipf_space(parser, "N", "the maximum space size of the Zipfian distribution which is used as the sampling method for the second node of one term in the path guided linear 1D SGD model (default: max path lengths)", {'k', "path-sgd-zipf-space"});
    args::ValueFlag<std::string> p_sgd_seed(parser, "STRING", "set the seed for the deterministic 1-threaded path guided linear 1D SGD model (default: pangenomic!)", {'q', "path-sgd-seed"});
    args::ValueFlag<std::string> p_sgd_snapshot(parser, "STRING", "set the prefix to which each snapshot graph of a path guided 1D SGD iteration should be written to, no default", {'u', "path-sgd-snapshot"});
    /// pipeline
    args::ValueFlag<std::string> pipeline(parser, "STRING", "apply a series of sorts, based on single-character command line arguments to this command, adding 's' as the default topological sort, 'f' to reverse the sort order, and 'g' to apply graph grooming", {'p', "pipeline"});
    /// paths
    args::Flag paths_by_min_node_id(parser, "paths-min", "sort paths by their lowest contained node id", {'L', "paths-min"});
    args::Flag paths_by_max_node_id(parser, "paths-max", "sort paths by their highest contained node id", {'M', "paths-max"});
    args::Flag paths_by_avg_node_id(parser, "paths-avg", "sort paths by their average contained node id", {'A', "paths-avg"});
    args::Flag paths_by_avg_node_id_rev(parser, "paths-avg-rev", "sort paths in reverse by their average contained node id", {'R', "paths-avg-rev"});
    args::ValueFlag<std::string> path_delim(parser, "path-delim", "sort paths in bins by their prefix up to this delimiter", {'D', "path-delim"});
    args::Flag progress(parser, "progress", "display progress of the sort", {'P', "progress"});
    args::Flag optimize(parser, "optimize", "use the MutableHandleGraph::optimize method", {'O', "optimize"});
    args::ValueFlag<uint64_t> nthreads(parser, "N", "number of threads to use for parallel sorters (currently only SGD is supported)", {'t', "threads"});

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
        std::cerr << "[odgi sort] Error: Please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (!dg_out_file) {
        std::cerr << "[odgi sort] Error: Please specify an output file to store the graph via -o=[FILE], --out=[FILE]." << std::endl;
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
    /*
    if (args::get(show_sort)) {
        std::vector<handle_t> order = (args::get(lazy) ? algorithms::lazy_topological_order(&graph) : algorithms::topological_order(&graph));
        for (auto& handle : order) {
            std::cout << graph.get_id(handle) << std::endl;
        }
    }
    */

    // default settings
    uint64_t df_chunk_size = args::get(depth_first_chunk) ? args::get(depth_first_chunk) : 1000;
    uint64_t bf_chunk_size = args::get(breadth_first_chunk) ? args::get(breadth_first_chunk) : std::numeric_limits<uint64_t>::max();
    uint64_t sgd_bandwidth = args::get(lsgd_bandwidth) ? args::get(lsgd_bandwidth) : 1000;
    double sgd_sampling_rate = args::get(lsgd_sampling_rate) ? args::get(lsgd_sampling_rate) : 20;
    double sgd_iter_max = args::get(lsgd_iter_max) ? args::get(lsgd_iter_max) : 30;
    double sgd_eps = args::get(lsgd_eps) ? args::get(lsgd_eps) : 0.01;
    double sgd_delta = args::get(lsgd_delta) ? args::get(lsgd_delta) : 0;
    uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;
    bool sgd_use_paths = args::get(lsgd_use_paths);
    /// path guided linear 1D SGD sort helpers
    // TODO beautify this, maybe put into its own file
    std::function<uint64_t(const std::vector<path_handle_t> &,
                           const xp::XP &)> get_sum_path_lengths
            = [&](const std::vector<path_handle_t> &path_sgd_use_paths, const xp::XP &path_index) {
                uint64_t sum_path_length = 0;
                for (auto& path : path_sgd_use_paths) {
                    sum_path_length += path_index.get_path_length(path);
                }
                return sum_path_length;
              };
    std::function<uint64_t(const std::vector<path_handle_t> &,
                           const xp::XP &)> get_max_path_length
            = [&](const std::vector<path_handle_t> &path_sgd_use_paths, const xp::XP &path_index) {
                uint64_t max_path_length = 0;
                for (auto& path : path_sgd_use_paths) {
                    max_path_length = std::max(max_path_length, path_index.get_path_length(path));
                }
                return max_path_length;
              };
    // default parameters
    std::string path_sgd_seed;
    if (p_sgd_seed) {
        if (num_threads > 1) {
            std::cerr << "[odgi sort] Error: Please only specify a seed for the path guided 1D linear SGD when using 1 thread." << std::endl;
            return 1;
        }
        path_sgd_seed = args::get(p_sgd_seed);
    } else {
        path_sgd_seed = "pangenomic!";
    }
    if (p_sgd_min_term_updates_paths && p_sgd_min_term_updates_num_nodes) {
        std::cerr << "[odgi sort] Error: There can only be on argument provided for the minimum number of term updates in the path guided 1D SGD."
                     "Please either use -G=[N], path-sgd-min-term-updates-paths=[N] or -U=[N], path-sgd-min-term-updates-nodes=[N]." << std::endl;
        return 1;
    }
    uint64_t path_sgd_iter_max = args::get(p_sgd_iter_max) ? args::get(p_sgd_iter_max) : 30;
    uint64_t path_sgd_iter_max_learning_rate = args::get(p_sgd_iter_with_max_learning_rate) ? args::get(p_sgd_iter_with_max_learning_rate) : 0;
    double path_sgd_zipf_theta = args::get(p_sgd_zipf_theta) ? args::get(p_sgd_zipf_theta) : 0.99;
    double path_sgd_eps = args::get(p_sgd_eps) ? args::get(p_sgd_eps) : 0.01;
    double path_sgd_delta = args::get(p_sgd_delta) ? args::get(p_sgd_delta) : 0;
    double path_sgd_max_eta = args::get(p_sgd_eta_max) ? args::get(p_sgd_eta_max) : graph.get_node_count();
    // will be filled, if the user decides to write a snapshot of the graph after each sorting iterationn
    std::vector<std::vector<handle_t>> snapshots;
    const bool snapshot = p_sgd_snapshot;
    const bool sample_from_nodes = p_sgd_sample_from_nodes;
    // default parameters that need a path index to be present
    uint64_t path_sgd_min_term_updates;
    uint64_t path_sgd_zipf_space;
    std::vector<path_handle_t> path_sgd_use_paths;
    xp::XP path_index;
    bool first_time_index = true;
    if (p_sgd || args::get(pipeline).find('Y') != std::string::npos) {
        // take care of path index
        if (xp_in_file) {
            std::ifstream in;
            in.open(args::get(xp_in_file));
            path_index.load(in);
            in.close();
        } else {
            path_index.from_handle_graph(graph);
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
                    std::cerr << "[odgi sort] Error: Path '" << buf
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
        uint64_t sum_path_length = get_sum_path_lengths(path_sgd_use_paths, path_index);
        if (args::get(p_sgd_min_term_updates_paths)) {
            path_sgd_min_term_updates = args::get(p_sgd_min_term_updates_paths) * sum_path_length;
        } else {
            if (args::get(p_sgd_min_term_updates_num_nodes)) {
                path_sgd_min_term_updates = args::get(p_sgd_min_term_updates_num_nodes) * graph.get_node_count();
            } else {
                path_sgd_min_term_updates = 0.1 * sum_path_length;
            }
        }
        path_sgd_zipf_space = args::get(p_sgd_zipf_space) ? args::get(p_sgd_zipf_space) : get_max_path_length(path_sgd_use_paths, path_index);
    }

    // helper, TODO: move into its own file
    // make a dagified copy, get its sort, and apply the order to our graph

    // did we groom the graph?
    bool was_groomed = false;
    std::string outfile = args::get(dg_out_file);
    if (outfile.size()) {
        if (args::get(eades)) {
            graph.apply_ordering(algorithms::eades_algorithm(&graph), true);
        } else if (args::get(two)) {
            graph.apply_ordering(algorithms::two_way_topological_order(&graph), true);
        } else if (args::get(optimize)) {
            graph.optimize();
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
        } else if (args::get(mondriaan)) {
            graph.apply_ordering(algorithms::mondriaan_sort(graph,
                                                            args::get(mondriaan_n_parts),
                                                            args::get(mondriaan_epsilon),
                                                            args::get(mondriaan_path_weight), false), true);
        } else if (args::get(lsgd)) {
            graph.apply_ordering(
                algorithms::linear_sgd_order(graph,
                                             sgd_bandwidth,
                                             sgd_sampling_rate,
                                             sgd_use_paths,
                                             sgd_iter_max,
                                             sgd_eps,
                                             sgd_delta,
                                             num_threads), true);
        } else if (args::get(p_sgd)) {
            std::vector<handle_t> order =
                algorithms::path_linear_sgd_order(graph,
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
                                                  num_threads,
                                                  progress,
                                                  path_sgd_seed,
                                                  snapshot,
                                                  snapshots,
                                                  sample_from_nodes);
            // TODO Check if we have to emit the snapshots
            if (snapshot) {
                std::string snapshot_prefix = args::get(p_sgd_snapshot);
                for (int j = 0; j < snapshots.size(); j++) {
                    std::cerr << "[path sgd sort]: Applying order to graph of iteration: " << std::to_string(j + 1) << std::endl;
                    std::string local_snapshot_prefix = snapshot_prefix + std::to_string(j + 1);
                    graph_t graph_copy = graph;
                    graph_copy.apply_ordering(snapshots[j], true);
                    ofstream f(local_snapshot_prefix);
                    graph_copy.serialize(f);
                    f.close();
                }
            }
            graph.apply_ordering(order, true);
        } else if (args::get(breadth_first)) {
            graph.apply_ordering(algorithms::breadth_first_topological_order(graph, bf_chunk_size), true);
        } else if (args::get(depth_first)) {
            graph.apply_ordering(algorithms::depth_first_topological_order(graph, df_chunk_size), true);
        } else if (args::get(randomize)) {
            graph.apply_ordering(algorithms::random_order(graph), true);
        } else if (!args::get(pipeline).empty()) {
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
                case 'e':
                    order = algorithms::eades_algorithm(&graph);
                    break;
                case 'd':
                {
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
                case 'S':
                    order = algorithms::linear_sgd_order(graph,
                                                         sgd_bandwidth,
                                                         sgd_sampling_rate,
                                                         sgd_use_paths,
                                                         sgd_iter_max,
                                                         sgd_eps,
                                                         sgd_delta,
                                                         num_threads);
                    break;
                case 'Y': {
                    if (!first_time_index) {
                        path_index.clean();
                        path_index.from_handle_graph(graph);
                    } else {
                        first_time_index = false;
                    }
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
                                                              num_threads,
                                                              progress,
                                                              path_sgd_seed,
                                                              snapshot,
                                                              snapshots,
                                                              sample_from_nodes);
                    break;
                }
                case 'f':
                    order.clear();
                    graph.for_each_handle([&order](const handle_t& handle) {
                            order.push_back(handle);
                        });
                    std::reverse(order.begin(), order.end());
                    break;
                case 'm':
                    order = algorithms::mondriaan_sort(graph,
                                                       args::get(mondriaan_n_parts),
                                                       args::get(mondriaan_epsilon),
                                                       args::get(mondriaan_path_weight), false);
                    break;
                case 'g': {
                    graph_t groomed;
                    algorithms::groom(graph, groomed);
                    graph = groomed;
                    was_groomed = true;
                    break;
                }
                default:
                    break;
                }
                if (!was_groomed) {
                    if (order.size() != graph.get_node_count()) {
                        std::cerr << "[odgi sort] Error: expected " << graph.get_node_count()
                                  << " handles in the order "
                                  << "but got " << order.size() << std::endl;
                        assert(false);
                    }
                    graph.apply_ordering(order, true);
                } else {
                    was_groomed = false;
                }
            }
        } else {
            graph.apply_ordering(algorithms::topological_order(&graph, true, false, args::get(progress)), true);
        }
        if (args::get(paths_by_min_node_id)) {
            graph.apply_path_ordering(algorithms::prefix_and_id_ordered_paths(graph, args::get(path_delim), false, false));
        }
        if (args::get(paths_by_max_node_id)) {
            graph.apply_path_ordering(algorithms::prefix_and_id_ordered_paths(graph, args::get(path_delim), false, true));
        }
        if (args::get(paths_by_avg_node_id)) {
            graph.apply_path_ordering(algorithms::prefix_and_id_ordered_paths(graph, args::get(path_delim), true, false));
        }
        if (args::get(paths_by_avg_node_id_rev)) {
            graph.apply_path_ordering(algorithms::prefix_and_id_ordered_paths(graph, args::get(path_delim), true, true));
        }
        if (outfile == "-") {
            graph.serialize(std::cout);
        } else {
            ofstream f(outfile.c_str());
            graph.serialize(f);
            f.close();
        }
    }
    return 0;
}

static Subcommand odgi_sort("sort", "sort a variation graph",
                              PIPELINE, 3, main_sort);


}
