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
        "establish 2D layouts of the graph using path-guided stochastic gradient descent (the graph must be sorted and id-compacted)");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> layout_out_file(parser, "FILE", "write the layout coordinates to this file in .lay binary format", {'o', "out"});
    args::ValueFlag<std::string> tsv_out_file(parser, "FILE", "write the TSV layout to this file", {'T', "tsv"});
    args::ValueFlag<std::string> svg_out_file(parser, "FILE", "write an SVG rendering to this file", {'s', "svg"});
    args::ValueFlag<std::string> png_out_file(parser, "FILE", "write a rasterized PNG rendering to this file", {'p', "png"});
    args::ValueFlag<uint64_t> png_height(parser, "FILE", "height of PNG rendering (default: 1000)", {'H', "png-height"});
    args::ValueFlag<uint64_t> png_border(parser, "FILE", "size of PNG border in bp (default: 10)", {'E', "png-border"});
    args::Flag color_paths(parser, "color-paths", "color paths (in PNG output)", {'C', "color-paths"});
    args::ValueFlag<double> render_scale(parser, "N", "image scaling (default 1.0)", {'R', "scale"});
    args::ValueFlag<double> render_border(parser, "N", "image border (in approximate bp) (default 100.0)", {'B', "border"});
    args::ValueFlag<std::string> xp_in_file(parser, "FILE", "load the path index from this file", {'X', "path-index"});
    /// Path-guided-2D-SGD parameters
    args::ValueFlag<std::string> p_sgd_in_file(parser, "FILE",
                                               "specify a line separated list of paths to sample from for the on the fly term generation process in the path guided linear 1D SGD default: sample from all paths)",
                                               {'f', "path-sgd-use-paths"});
    args::ValueFlag<char> p_sgd_layout_initialization(parser, "C", "specify the layout initialization mode: d) node rank in X and gaussian noise in Y (default)\nu) node rank in X and uniform noise in Y\ng) gaussian noise in X and Y\nh) hilbert curve in X and Y", {'N', "layout-initialization"});
    args::ValueFlag<double> p_sgd_min_term_updates_paths(parser, "N",
                                                         "minimum number of terms to be updated before a new path guided linear 1D SGD iteration with adjusted learning rate eta starts, expressed as a multiple of total path length (default: 10)",
                                                         {'G', "path-sgd-min-term-updates-paths"});
    args::ValueFlag<double> p_sgd_min_term_updates_num_nodes(parser, "N",
                                                             "minimum number of terms to be updated before a new path guided linear 1D SGD iteration with adjusted learning rate eta starts, expressed as a multiple of the number of nodes (default: argument is not set, the default of -G=[N], path-sgd-min-term-updates-paths=[N] is used)",
                                                             {'U', "path-sgd-min-term-updates-nodes"});
    args::ValueFlag<double> p_sgd_delta(parser, "N",
                                        "threshold of maximum displacement approximately in bp at which to stop path guided linear 1D SGD (default: 0)",
                                        {'j', "path-sgd-delta"});
    args::ValueFlag<double> p_sgd_eps(parser, "N",
                                      "final learning rate for path guided linear 1D SGD model (default: 0.01)",
                                      {'g', "path-sgd-eps"});
    args::ValueFlag<double> p_sgd_eta_max(parser, "N",
                                          "first and maximum learning rate for path guided linear 1D SGD model (default: squared longest path length)",
                                          {'v', "path-sgd-eta-max"});
    args::ValueFlag<double> p_sgd_zipf_theta(parser, "N",
                                             "the theta value for the Zipfian distribution which is used as the sampling method for the second node of one term in the path guided linear 1D SGD model (default: 0.99)",
                                             {'a', "path-sgd-zipf-theta"});
    args::ValueFlag<uint64_t> p_sgd_iter_max(parser, "N",
                                             "max number of iterations for path guided linear 1D SGD model (default: 30)",
                                             {'x', "path-sgd-iter-max"});
    args::ValueFlag<uint64_t> p_sgd_iter_with_max_learning_rate(parser, "N",
                                                                "iteration where the learning rate is max for path guided linear 1D SGD model (default: 0)",
                                                                {'F', "iteration-max-learning-rate"});
    args::ValueFlag<uint64_t> p_sgd_zipf_space(parser, "N",
                                               "the maximum space size of the Zipfian distribution which is used as the sampling method for the second node of one term in the path guided linear 1D SGD model (default: min(max path lengths, 10000))",
                                               {'k', "path-sgd-zipf-space"});
    args::ValueFlag<uint64_t> p_sgd_zipf_space_max(parser, "N", "the maximum space size of the Zipfian distribution beyond which quantization occurs (default: 1000)", {'I', "path-sgd-zipf-space-max"});
    args::ValueFlag<uint64_t> p_sgd_zipf_space_quantization_step(parser, "N", "quantization step when the maximum space size of the Zipfian distribution is exceeded (default: 100)", {'l', "path-sgd-zipf-space-quantization-step"});

    args::ValueFlag<std::string> p_sgd_seed(parser, "STRING",
                                            "set the seed for the deterministic 1-threaded path guided linear 1D SGD model (default: pangenomic!)",
                                            {'q', "path-sgd-seed"});
    args::ValueFlag<std::string> p_sgd_snapshot(parser, "STRING",
                                                "set the prefix to which each snapshot graph of a path guided 1D SGD iteration should be written to, no default",
                                                {'u', "path-sgd-snapshot"});

    //args::ValueFlag<double> x_pad(parser, "N", "padding between connected component layouts (default 10.0)",
    //{'p', "x-padding"});
    args::Flag progress(parser, "progress", "display progress of the sort", {'P', "progress"});
    args::ValueFlag<uint64_t> nthreads(parser, "N",
                                       "number of threads to use for parallel sorters (currently only SGD is supported)",
                                       {'t', "threads"});
    args::Flag debug(parser, "debug", "print information about the layout", {'d', "debug"});

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
            << "[odgi::layout] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
            << std::endl;
        return 1;
    }

    if (!layout_out_file && !svg_out_file && !png_out_file && !tsv_out_file) {
        std::cerr
            << "[odgi::layout] error: please specify an output file to where to store the layout via -o/--out=[FILE], -p/--png=[FILE], -s/--svg=[FILE], -T/--tsv=[FILE]"
            << std::endl;
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

    const uint64_t t_max = !p_sgd_iter_max ? 30 : args::get(p_sgd_iter_max);
    const double eps = !p_sgd_eps ? 0.01 : args::get(p_sgd_eps);
    //const double x_padding = !x_pad ? 10.0 : args::get(x_pad);
    const double svg_scale = !render_scale ? 1.0 : args::get(render_scale);
    const double border_bp = !render_border ? 100.0 : args::get(render_border);
    const double sgd_delta = p_sgd_delta ? args::get(p_sgd_delta) : 0;
    const uint64_t num_threads = nthreads ? args::get(nthreads) : 1;
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
              uint64_t max_path_step_count = 0;
              for (auto& path : path_sgd_use_paths) {
                  max_path_step_count = std::max(max_path_step_count, path_index.get_path_step_count(path));
              }
              return max_path_step_count;
          };
    // default parameters
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
    path_sgd_zipf_space = args::get(p_sgd_zipf_space) ? std::min(args::get(p_sgd_zipf_space), max_path_step_count) : std::min((uint64_t)10000, max_path_step_count);
    double path_sgd_max_eta = args::get(p_sgd_eta_max) ? args::get(p_sgd_eta_max) : (double) max_path_step_count * max_path_step_count;

    path_sgd_zipf_space_max = args::get(p_sgd_zipf_space_max) ? std::min(path_sgd_zipf_space, args::get(p_sgd_zipf_space_max)) : 1000;
    path_sgd_zipf_space_quantization_step = args::get(p_sgd_zipf_space_quantization_step) ? std::max((uint64_t)2, args::get(p_sgd_zipf_space_quantization_step)) : 100;

    /*// refine order by weakly connected components
      std::vector<ska::flat_hash_set<handlegraph::nid_t>> weak_components = algorithms::weakly_connected_components(
      &graph);
      #ifdef debug_components
      std::cerr << "components count: " << weak_components.size() << std::endl;
      #endif
      std::vector<std::pair<double, uint64_t>> weak_component_order;
      for (int i = 0; i < weak_components.size(); i++) {
      auto &weak_component = weak_components[i];
      uint64_t id_sum = 0;
      for (auto node_id : weak_component) {
      id_sum += node_id;
      }
      double avg_id = id_sum / (double) weak_component.size();
      weak_component_order.push_back(std::make_pair(avg_id, i));
      }
      std::sort(weak_component_order.begin(), weak_component_order.end());*/

    std::vector<std::atomic<double>> graph_X(graph.get_node_count() * 2);  // Graph's X coordinates for node+ and node-
    std::vector<std::atomic<double>> graph_Y(graph.get_node_count() * 2);  // Graph's Y coordinates for node+ and node-

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> uniform_noise(0, sqrt(graph.get_node_count() * 2));
    std::normal_distribution<double> gaussian_noise(0,  sqrt(graph.get_node_count() * 2));

    uint64_t len = 0;
    nid_t last_node_id = graph.min_node_id();
    char layout_initialization = p_sgd_layout_initialization ? args::get(p_sgd_layout_initialization) : 'd';

    uint64_t square_space = graph.get_node_count() * 2;
    uint64_t x, y;
    graph.for_each_handle([&](const handle_t &h) {
          nid_t node_id = graph.get_id(h);
          if (node_id - last_node_id > 1) {
              std::cerr << "[odgi::layout] error: the graph is not optimized. Please run 'odgi sort' using -O, --optimize" << std::endl;
              exit(1);
          }
          last_node_id = node_id;

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
        num_threads,
        show_progress,
        snapshot,
        snapshot_prefix,
        graph_X,
        graph_Y
        );

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

    if (tsv_out_file) {
        auto& outfile = args::get(tsv_out_file);
        if (outfile.size()) {
            if (outfile == "-") {
                algorithms::layout::to_tsv(std::cout, X_final, Y_final, graph);
            } else {
                ofstream f(outfile.c_str());
                algorithms::layout::to_tsv(f, X_final, Y_final, graph);
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

    if (svg_out_file) {
        auto& outfile = args::get(svg_out_file);
        ofstream f(outfile.c_str());
        algorithms::draw_svg(f, X_final, Y_final, graph, svg_scale, border_bp);
        f.close();    
    }

    if (png_out_file) {
        auto& outfile = args::get(png_out_file);
        uint64_t _png_height = png_height ? args::get(png_height) : 1000;
        bool _color_paths = args::get(color_paths);
        algorithms::draw_png(outfile, X_final, Y_final, graph, 1.0, border_bp, 0, _png_height, 0.0, 1.0, _color_paths);
    }
    
    return 0;
}

static Subcommand odgi_layout("layout", "use path guided SGD to make 2D layouts of the graph",
                               PIPELINE, 3, main_layout);


}
