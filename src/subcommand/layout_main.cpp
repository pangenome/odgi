#include "subcommand.hpp"
#include <iostream>
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/xp.hpp"
#include "algorithms/sgd_layout.hpp"
#include "algorithms/path_sgd_layout.hpp"
#include "algorithms/draw.hpp"
#include "algorithms/layout.hpp"
#include "hilbert.hpp"
#include "utils.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_layout(int argc, char **argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc - 1; ++i) {
        argv[i] = argv[i + 1];
    }
    std::string prog_name = "odgi layout";
    argv[0] = (char *) prog_name.c_str();
    --argc;

    args::ArgumentParser parser(
        "Establish 2D layouts of the graph using path-guided stochastic gradient descent. The graph must be sorted and id-compacted.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::Group files_io_opts(parser, "[ Files IO ]");
    args::ValueFlag<std::string> layout_out_file(files_io_opts, "FILE", "Write the layout coordinates to this FILE in .lay binary format.", {'o', "out"});
    args::ValueFlag<std::string> tsv_out_file(files_io_opts, "FILE", "Write the layout in TSV format to this FILE.", {'T', "tsv"});
    args::ValueFlag<std::string> xp_in_file(files_io_opts, "FILE", "Load the path index from this FILE so that it does not have to be created for the layout calculation.", {'X', "path-index"});
    args::ValueFlag<std::string> tmp_base(files_io_opts, "PATH", "directory for temporary files", {'C', "temp-dir"});
    /// Path-guided-2D-SGD parameters
    args::ValueFlag<std::string> p_sgd_in_file(files_io_opts, "FILE",
                                               "Specify a line separated list of paths to sample from for the on the fly term generation process in the path guided 2D SGD (default: sample from all paths).",
                                               {'f', "path-sgd-use-paths"});
    args::Group layout_init_opts(parser, "[ Layout Initialization Options ]");
    args::ValueFlag<char> p_sgd_layout_initialization(layout_init_opts, "C", "Specify the layout initialization mode:\nd) Node rank in X and gaussian noise in Y (default).\nr) Uniform noise in X and Y in the order of the graph length.\nu) Node rank in X and uniform noise in Y.\ng) Gaussian noise in X and Y.\nh) Hilbert curve in X and Y.", {'N', "layout-initialization"});
    args::Group pg_sgd_opts(parser, "[ PG-SGD Options ]");
    args::ValueFlag<double> p_sgd_min_term_updates_paths(pg_sgd_opts, "N",
                                                         "Minimum number of terms N to be updated before a new path guided 2D SGD iteration with adjusted learning rate eta starts, expressed as a multiple of total path length (default: 10).",
                                                         {'G', "path-sgd-min-term-updates-paths"});
    args::ValueFlag<double> p_sgd_min_term_updates_num_nodes(pg_sgd_opts, "N",
                                                             "Minimum number of terms N to be updated before a new path guided linear 1D SGD iteration with adjusted learning rate eta starts, expressed as a multiple of the number of nodes (default: argument is not set, the default of -G=[N], path-sgd-min-term-updates-paths=[N] is used).",
                                                             {'U', "path-sgd-min-term-updates-nodes"});
    args::ValueFlag<double> p_sgd_delta(pg_sgd_opts, "N",
                                        "The threshold of the maximum displacement approximately in bp at which to stop path guided 2D SGD (default: 0).",
                                        {'j', "path-sgd-delta"});
    args::ValueFlag<double> p_sgd_eps(pg_sgd_opts, "N",
                                      "The final learning rate for path guided 2D SGD model (default: 0.01).",
                                      {'g', "path-sgd-eta"});
    args::ValueFlag<double> p_sgd_eta_max(pg_sgd_opts, "N",
                                          "The first and maximum learning rate N for path guided 2D SGD model (default: squared longest path length).",
                                          {'v', "path-sgd-eta-max"});
    args::ValueFlag<double> p_sgd_zipf_theta(pg_sgd_opts, "N",
                                             "The theta value N for the Zipfian distribution which is used as the sampling method for the second node of one term in the path guided 2D SGD model (default: 0.99).",
                                             {'a', "path-sgd-zipf-theta"});
    args::ValueFlag<uint64_t> p_sgd_iter_max(pg_sgd_opts, "N",
                                             "The maximum number of iterations N for the path guided 2D SGD model (default: 30).",
                                             {'x', "path-sgd-iter-max"});
    args::ValueFlag<double> p_sgd_cooling(pg_sgd_opts, "N",
                                          "Use this fraction of the iterations for layout annealing (default: 0.5).",
                                          {'K', "path-sgd-cooling"});
    args::ValueFlag<uint64_t> p_sgd_iter_with_max_learning_rate(pg_sgd_opts, "N",
                                                                "Specify the iteration N where the learning rate is max for path guided 2D SGD model (default: 0).",
                                                                {'F', "path-sgd-iteration-max-learning-rate"});
    args::ValueFlag<uint64_t> p_sgd_zipf_space(pg_sgd_opts, "N",
                                               "The maximum space size N of the Zipfian distribution which is used as the sampling method for the second node of one term in the path guided 2D SGD model (default: max path lengths).",
                                               {'k', "path-sgd-zipf-space"});
    args::ValueFlag<uint64_t> p_sgd_zipf_space_max(pg_sgd_opts, "N", "The maximum space size N of the Zipfian distribution beyond which quantization occurs (default: 1000).", {'I', "path-sgd-zipf-space-max"});
    args::ValueFlag<uint64_t> p_sgd_zipf_space_quantization_step(pg_sgd_opts, "N", "The size of the quantization step N when the maximum space size of the Zipfian distribution is exceeded (default: 100).", {'l', "path-sgd-zipf-space-quantization-step"});
    /*
    args::ValueFlag<std::string> p_sgd_seed(parser, "STRING",
                                            "set the seed for the deterministic 1-threaded path guided linear 1D SGD model (default: pangenomic!)",
                                            {'q', "path-sgd-seed"});
    */
    args::ValueFlag<std::string> p_sgd_snapshot(pg_sgd_opts, "STRING",
                                                "Set the prefix to which each snapshot layout of a path guided 2D SGD iteration should be written to (default: NONE).",
                                                {'u', "path-sgd-snapshot"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> nthreads(threading_opts, "N",
                                       "Number of threads to use for parallel operations.",
                                       {'t', "threads"});
#ifdef USE_GPU
    // GPU-enabled Layout
    args::Group gpu_opts(parser, "[ GPU ]");
    args::Flag gpu_compute(gpu_opts, "gpu", "Enable computation with GPU.", {"gpu"});
#endif
    args::Group processing_info_opts(parser, "[ Processsing Information ]");
    args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help summary for odgi layout.", {'h', "help"});

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

    if (!dg_in_file) {
        std::cerr
            << "[odgi::layout] error: Please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
            << std::endl;
        return 1;
    }

    if (!layout_out_file && !tsv_out_file) {
        std::cerr
            << "[odgi::layout] error: Please specify an output file to where to store the layout via -o/--out=[FILE] or -T/--tsv=[FILE]."
            << std::endl;
        return 1;
    }

	const uint64_t num_threads = nthreads ? args::get(nthreads) : 1;

	graph_t graph;
    assert(argc > 0);
    if (!args::get(dg_in_file).empty()) {
        std::string infile = args::get(dg_in_file);
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
			utils::handle_gfa_odgi_input(infile, "layout", args::get(progress), num_threads, graph);
        }
    }

    if (tmp_base) {
        xp::temp_file::set_dir(args::get(tmp_base));
    } else {
        char cwd[512];
        getcwd(cwd, sizeof(cwd));
        xp::temp_file::set_dir(std::string(cwd));
    }

    if (!graph.is_optimized()) {
		std::cerr << "[odgi::layout] error: the graph is not optimized. Please run 'odgi sort' using -O, --optimize." << std::endl;
		exit(1);
    }

    const double eps = !p_sgd_eps ? 0.01 : args::get(p_sgd_eps);
    const double sgd_delta = p_sgd_delta ? args::get(p_sgd_delta) : 0;
    const bool show_progress = progress ? args::get(progress) : false;
    /// path guided linear 2D SGD sort helpers
    // TODO beautify this, maybe put into its own file
    std::function<uint64_t(const std::vector<path_handle_t> &,
                           const xp::XP &)> get_sum_path_step_count
        = [&](const std::vector<path_handle_t> &path_sgd_use_paths, const xp::XP &path_index) {
              uint64_t sum_path_step_count = 0;
              for (auto& path : path_sgd_use_paths) {
                  sum_path_step_count += path_index.get_path_step_count(path);
              }
              return sum_path_step_count;
          };
    std::function<uint64_t(const std::vector<path_handle_t> &,
                           const xp::XP &)> get_max_path_step_count
        = [&](const std::vector<path_handle_t> &path_sgd_use_paths, const xp::XP &path_index) {
              size_t max_path_step_count = 0;
              for (auto& path : path_sgd_use_paths) {
                  max_path_step_count = std::max(max_path_step_count, path_index.get_path_step_count(path));
              }
              return max_path_step_count;
          };
    // default parameters
    /* We don't do this, yet.
    std::string path_sgd_seed;
    if (p_sgd_seed) {
        if (num_threads > 1) {
            std::cerr
                << "[odgi::layout] error: please only specify a seed for the path guided 2D linear SGD when using 1 thread."
                << std::endl;
            return 1;
        }
        path_sgd_seed = args::get(p_sgd_seed);
    } else {
        path_sgd_seed = "pangenomic!";
    }
    */
    if (p_sgd_min_term_updates_paths && p_sgd_min_term_updates_num_nodes) {
        std::cerr
            << "[odgi::layout] error: there can only be one argument provided for the minimum number of term updates in the path guided 1D SGD."
            "Please either use -G=[N], path-sgd-min-term-updates-paths=[N] or -U=[N], path-sgd-min-term-updates-nodes=[N]."
            << std::endl;
        return 1;
    }
    uint64_t path_sgd_iter_max = p_sgd_iter_max ? args::get(p_sgd_iter_max) : 30;
    uint64_t path_sgd_iter_max_learning_rate = p_sgd_iter_with_max_learning_rate ? args::get(
        p_sgd_iter_with_max_learning_rate) : 0;
    double path_sgd_zipf_theta = p_sgd_zipf_theta ? args::get(p_sgd_zipf_theta) : 0.99;
    double path_sgd_eps = p_sgd_eps ? args::get(p_sgd_eps) : 0.01;
    double path_sgd_delta = p_sgd_delta ? args::get(p_sgd_delta) : 0;
    double path_sgd_cooling = p_sgd_cooling ? args::get(p_sgd_cooling) : 0.5;
    // will be filled, if the user decides to write a snapshot of the graph after each sorting iterationn
    const bool snapshot = p_sgd_snapshot;
    std::string snapshot_prefix;
    if (snapshot) {
        snapshot_prefix = args::get(p_sgd_snapshot);
    }

    // default parameters that need a path index to be present
    uint64_t path_sgd_min_term_updates;
    uint64_t path_sgd_zipf_space, path_sgd_zipf_space_max, path_sgd_zipf_space_quantization_step;
    std::vector<path_handle_t> path_sgd_use_paths;
    xp::XP path_index;
    bool first_time_index = true;

    // take care of path index
    if (xp_in_file) {
        std::ifstream in;
        in.open(args::get(xp_in_file));
        path_index.load(in);
        in.close();
    } else {
        path_index.from_handle_graph(graph, num_threads);
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
                std::cerr << "[odgi::layout] error: path '" << buf
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


    uint64_t sum_path_step_count = get_sum_path_step_count(path_sgd_use_paths, path_index);
    if (args::get(p_sgd_min_term_updates_paths)) {
        path_sgd_min_term_updates = args::get(p_sgd_min_term_updates_paths) * sum_path_step_count;
    } else {
        if (args::get(p_sgd_min_term_updates_num_nodes)) {
            path_sgd_min_term_updates = args::get(p_sgd_min_term_updates_num_nodes) * graph.get_node_count();
        } else {
            path_sgd_min_term_updates = 10.0 * sum_path_step_count;
        }
    }
    uint64_t max_path_step_count = get_max_path_step_count(path_sgd_use_paths, path_index);
    path_sgd_zipf_space = args::get(p_sgd_zipf_space) ? std::min(args::get(p_sgd_zipf_space), max_path_step_count) : max_path_step_count;
    double path_sgd_max_eta = args::get(p_sgd_eta_max) ? args::get(p_sgd_eta_max) : (double) max_path_step_count * max_path_step_count;

    path_sgd_zipf_space_max = args::get(p_sgd_zipf_space_max) ? std::min(path_sgd_zipf_space, args::get(p_sgd_zipf_space_max)) : 1000;
    path_sgd_zipf_space_quantization_step = args::get(p_sgd_zipf_space_quantization_step) ? std::max((uint64_t)2, args::get(p_sgd_zipf_space_quantization_step)) : 100;

    std::vector<std::atomic<double>> graph_X(graph.get_node_count() * 2);  // Graph's X coordinates for node+ and node-
    std::vector<std::atomic<double>> graph_Y(graph.get_node_count() * 2);  // Graph's Y coordinates for node+ and node-

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> uniform_noise(0, sqrt(graph.get_node_count() * 2));
    std::normal_distribution<double> gaussian_noise(0,  sqrt(graph.get_node_count() * 2));
    uint64_t total_length = 0;
    graph.for_each_handle([&](const handle_t &h) {
                              total_length += graph.get_length(h);
                          });
    std::uniform_real_distribution<double> uniform_noise_in_length(0, total_length);

    uint64_t len = 0;
    nid_t last_node_id = graph.min_node_id();
    char layout_initialization = p_sgd_layout_initialization ? args::get(p_sgd_layout_initialization) : 'd';

    uint64_t square_space = graph.get_node_count() * 2;
    uint64_t x, y;
    graph.for_each_handle([&](const handle_t &h) {
          uint64_t pos = 2 * number_bool_packing::unpack_number(h);
          switch (layout_initialization) {
          case 'g': {
              graph_X[pos].store(gaussian_noise(rng));
              graph_Y[pos].store(gaussian_noise(rng));
              graph_X[pos + 1].store(gaussian_noise(rng));
              graph_Y[pos + 1].store(gaussian_noise(rng));
              break;
          }
          case 'u': {
              graph_X[pos].store(len);
              graph_Y[pos].store(uniform_noise(rng));
              len += graph.get_length(h);
              graph_X[pos + 1].store(len);
              graph_Y[pos + 1].store(uniform_noise(rng));
              break;
          }
          case 'r': {
              graph_X[pos].store(uniform_noise_in_length(rng));
              graph_Y[pos].store(uniform_noise_in_length(rng));
              graph_X[pos + 1].store(uniform_noise_in_length(rng));
              graph_Y[pos + 1].store(uniform_noise_in_length(rng));
              break;
          }
          case 'h': {
              d2xy(square_space, pos, &x, &y);
              graph_X[pos].store(x);
              graph_Y[pos].store(y);
              d2xy(square_space, pos + 1, &x, &y);
              graph_X[pos + 1].store(x);
              graph_Y[pos + 1].store(y);
              break;
          }
          default: {
              graph_X[pos].store(len);
              graph_Y[pos].store(gaussian_noise(rng));
              len += graph.get_length(h);
              graph_X[pos + 1].store(len);
              graph_Y[pos + 1].store(gaussian_noise(rng));
          }
          }
          //std::cerr << pos << ": " << graph_X[pos] << "," << graph_Y[pos] << " ------ " << graph_X[pos + 1] << "," << graph_Y[pos + 1] << std::endl;
      });

    //double max_x = 0;
#ifdef USE_GPU
    if (gpu_compute) { // run on GPU
        algorithms::path_linear_sgd_layout_gpu(
            graph,
            path_index,
            path_sgd_use_paths,
            path_sgd_iter_max,
            0,
            path_sgd_min_term_updates,
            sgd_delta,
            eps,
            path_sgd_max_eta,
            path_sgd_zipf_theta,
            path_sgd_zipf_space,
            path_sgd_zipf_space_max,
            path_sgd_zipf_space_quantization_step,
            path_sgd_cooling,
            num_threads,
            show_progress,
            snapshot,
            snapshot_prefix,
            graph_X,
            graph_Y
            );
    } 
#endif

#ifdef USE_GPU
    if (!gpu_compute) { // run on CPU
#endif
        algorithms::path_linear_sgd_layout(
            graph,
            path_index,
            path_sgd_use_paths,
            path_sgd_iter_max,
            0,
            path_sgd_min_term_updates,
            sgd_delta,
            eps,
            path_sgd_max_eta,
            path_sgd_zipf_theta,
            path_sgd_zipf_space,
            path_sgd_zipf_space_max,
            path_sgd_zipf_space_quantization_step,
            path_sgd_cooling,
            num_threads,
            show_progress,
            snapshot,
            snapshot_prefix,
            graph_X,
            graph_Y
            );
#ifdef USE_GPU
    }
#endif
    // drop out of atomic stuff... maybe not the best way to do this
    // TODO: use directly the atomic vector?
    std::vector<double> X_final(graph_X.size());
    uint64_t i = 0;
    for (auto& x : graph_X) {
        X_final[i++] = x.load();
    }
    std::vector<double> Y_final(graph_Y.size());
    i = 0;
    for (auto& y : graph_Y) {
        Y_final[i++] = y.load();
    }

    // refine order by weakly connected components
    std::vector<std::vector<handlegraph::handle_t>> weak_components = algorithms::weakly_connected_component_vectors(&graph);

    //uint64_t num_components_on_each_dimension = std::ceil(sqrt(weak_components.size()));
    //std::cerr << " num_components_on_each_dimension " << num_components_on_each_dimension << std::endl;

    double border = 1000.0;
    double curr_y_offset = border;
    std::vector<algorithms::coord_range_2d_t> component_ranges;
    for (auto& component : weak_components) {
        component_ranges.emplace_back();
        auto& component_range = component_ranges.back();
        for (auto& handle : component) {
            uint64_t pos = 2 * number_bool_packing::unpack_number(handle);
            for (uint64_t j = pos; j <= pos+1; ++j) {
                component_range.include(X_final[j], Y_final[j]);
            }
        }
        component_range.x_offset = component_range.min_x - border;
        component_range.y_offset = curr_y_offset - component_range.min_y;
        curr_y_offset += component_range.height() + border;
    }

    for (uint64_t num_component = 0; num_component < weak_components.size(); ++num_component) {
        auto& component_range = component_ranges[num_component];

        for (auto& handle :  weak_components[num_component]) {
            uint64_t pos = 2 * number_bool_packing::unpack_number(handle);

            for (uint64_t j = pos; j <= pos+1; ++j) {
                X_final[j] -= component_range.x_offset;
                Y_final[j] += component_range.y_offset;
            }
        }
    }


    if (tsv_out_file) {
        auto& outfile = args::get(tsv_out_file);
        if (outfile.size()) {
            if (outfile == "-") {
                algorithms::layout::to_tsv(std::cout, X_final, Y_final, weak_components);
            } else {
                ofstream f(outfile.c_str());
                algorithms::layout::to_tsv(f, X_final, Y_final, weak_components);
                f.close();
            }
        }
    }

    if (layout_out_file) {
        auto& outfile = args::get(layout_out_file);
        if (outfile.size()) {
            algorithms::layout::Layout lay(X_final, Y_final);
            if (outfile == "-") {
                lay.serialize(std::cout);
            } else {
                ofstream f(outfile.c_str());
                lay.serialize(f);
                f.close();
            }
        }
    }
    
    return 0;
}

static Subcommand odgi_layout("layout", "Establish 2D layouts of the graph using path-guided stochastic gradient descent.",
                               PIPELINE, 3, main_layout);


}
