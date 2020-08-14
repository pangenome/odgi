#include "subcommand.hpp"
#include <iostream>
#include "odgi.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "algorithms/xp.hpp"
#include "algorithms/sgd_layout.hpp"
#include "algorithms/path_sgd_layout.hpp"

namespace odgi {

    using namespace odgi::subcommand;


    void draw_svg(ostream &out, const std::vector<double> &X, const std::vector<double> &Y, const HandleGraph &graph, double scale) {
        double border = 10.0;
        double min_x = 0;
        double min_y = 0;
        double max_x = 0;
        double max_y = 0;
        // determine boundaries
        uint64_t n = graph.get_node_count() * 2;
        for (uint64_t i = 0; i < n; ++i) {
            double x = X[i] * scale;
            double y = Y[i] * scale;
            if (x < min_x) min_x = x;
            if (x > max_x) max_x = x;
            if (y < min_y) min_y = y;
            if (y > max_y) max_y = y;
            std::cerr << "X" << i << ": " << X[i] << "Y" << i << ": " << Y[i] << std::endl;
        }

        double width = max_x - min_x;
        double height = max_y - min_y;
        std::cerr << "width: " << width << std::endl;
        std::cerr << "height: " << height << std::endl;

        out << "<svg width=\"" << width + border << "\" height=\"" << height + border << "\" "
            << "viewBox=\"" << min_x - border / 2 << " " << min_y - border / 2 << " " << width + border << " "
            << height + border << "\" xmlns=\"http://www.w3.org/2000/svg\">"
            << "<style type=\"text/css\">"
            << "line{stroke:black;stroke-width:1.0;stroke-opacity:1.0;stroke-linecap:round;}"
            //<< "circle{{r:" << 1.0 << ";fill:black;fill-opacity:" << 1.0 << "}}"
            << "</style>"
            << std::endl;

        graph.for_each_handle([&](const handle_t& handle) {
            uint64_t a = 2 * number_bool_packing::unpack_number(handle);

            //std::cerr << a << ": " << X[a] << "," << Y[a] << " ------ " << X[a + 1] << "," << Y[a + 1] << std::endl;
            out << "<line x1=\"" << X[a] * scale << "\" x2=\"" << X[a + 1] * scale
                << "\" y1=\"" << Y[a] * scale << "\" y2=\"" << Y[a + 1] * scale << "\"/>"
                << std::endl;
        });

        // iterate through graph edges
        /*graph.for_each_edge([&](const edge_t &e) {
            uint64_t a = 2 * number_bool_packing::unpack_number(e.first);
            uint64_t b = 2 * number_bool_packing::unpack_number(e.second);

            //std::cerr << a << ": " << X[a] << "," << Y[a] << " ------ " << X[a + 1] << "," << Y[a + 1] << std::endl;
            out << "<line x1=\"" << X[a] * scale << "\" x2=\"" << X[a + 1] * scale
                << "\" y1=\"" << Y[a] * scale << "\" y2=\"" << Y[a + 1] * scale << "\"/>"
                << std::endl;

            //std::cerr << b << ": " << X[b] << "," << Y[b] << " ------ " << X[b + 1] << "," << Y[b + 1] << std::endl;
            out << "<line x1=\"" << X[b] * scale << "\" x2=\"" << X[b + 1] * scale
                << "\" y1=\"" << Y[b] * scale << "\" y2=\"" << Y[b + 1] * scale << "\"/>"
                << std::endl;
        });*/


        /* // to draw nodes
        for (uint64_t i = 0; i < n; ++i) {
            std::cout << "<circle cx=\"" << X[i*2]*scale << "\" cy=\"" << X[i*2+1]*scale << "\" r=\"1.0\"/>" << std::endl;
        }
        */

        out << "</svg>" << std::endl;
    }

        void to_tsv(ostream &out, const std::vector<double> &X, const std::vector<double> &Y, const HandleGraph &graph, double scale) {
        // determine boundaries
        uint64_t n = graph.get_node_count() * 2;
        out << "idx" << "\t" << "X" << "\t" << "Y" << std::endl;
        for (uint64_t i = 0; i < n; ++i) {
            out << i << "\t" << X[i] << "\t" << Y[i] << std::endl;
            // std::cerr << "X" << i << ": " << X[i] << "Y" << i << ": " << Y[i] << std::endl;
        }
    }

    int main_layout2(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi layout2";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser(
                "draw 2D layouts of the graph using stocastic gradient descent (the graph must be sorted and id-compacted)");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
        args::ValueFlag<std::string> svg_out_file(parser, "FILE", "write the SVG rendering to this file",
                                                  {'o', "out"});
        args::ValueFlag<std::string> xp_in_file(parser, "FILE", "load the path index from this file",
                                                {'X', "path-index"});
        /// Path-guided-2D-SGD parameters
        args::Flag p_sgd_deterministic(parser, "path-sgd-deterministic",
                                       "run the path guided 1D linear SGD in deterministic mode, will automatically set the number of threads to 1, multithreading is not supported in this mode (default: flag not set)",
                                       {'I', "path-sgd-deterministic"});
        args::ValueFlag<std::string> p_sgd_in_file(parser, "FILE",
                                                   "specify a line separated list of paths to sample from for the on the fly term generation process in the path guided linear 1D SGD (default: sample from all paths)",
                                                   {'f', "path-sgd-use-paths"});
        args::ValueFlag<double> p_sgd_min_term_updates_paths(parser, "N",
                                                             "minimum number of terms to be updated before a new path guided linear 1D SGD iteration with adjusted learning rate eta starts, expressed as a multiple of total path length (default: 0.1)",
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
                                              "first and maximum learning rate for path guided linear 1D SGD model (default: number of nodes in the graph)",
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
                                                   "the maximum space size of the Zipfian distribution which is used as the sampling method for the second node of one term in the path guided linear 1D SGD model (default: max path lengths)",
                                                   {'k', "path-sgd-zipf-space"});
        args::ValueFlag<std::string> p_sgd_seed(parser, "STRING",
                                                "set the seed for the deterministic 1-threaded path guided linear 1D SGD model (default: pangenomic!)",
                                                {'q', "path-sgd-seed"});
        args::ValueFlag<std::string> p_sgd_snapshot(parser, "STRING",
                                                    "set the prefix to which each snapshot graph of a path guided 1D SGD iteration should be written to, no default",
                                                    {'u', "path-sgd-snapshot"});

        args::ValueFlag<double> x_pad(parser, "N", "padding between connected component layouts (default 10.0)",
                                      {'p', "x-padding"});
        args::ValueFlag<double> render_scale(parser, "N", "image scaling (default 5.0)", {'R', "render-scale"});
        args::ValueFlag<double> render_png(parser, "N", "PNG", {'n', "render-png"});

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
                    << "[odgi layout2] error: Please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
                    << std::endl;
            return 1;
        }

        if (!svg_out_file) {
            std::cerr
                    << "[odgi layout2] error: Please specify an output file to where to store the layout via -o=[FILE], --out=[FILE]."
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

        uint64_t t_max = !args::get(p_sgd_iter_max) ? 30 : args::get(p_sgd_iter_max);

        double eps = !args::get(p_sgd_eps) ? 0.01 : args::get(p_sgd_eps);
        double x_padding = !args::get(x_pad) ? 10.0 : args::get(x_pad);
        double svg_scale = !args::get(render_scale) ? 5.0 : args::get(render_scale);
        double sgd_delta = args::get(p_sgd_delta) ? args::get(p_sgd_delta) : 0;
        uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;
        /// path guided linear 2D SGD sort helpers
        // TODO beautify this, maybe put into its own file
        std::function<uint64_t(const std::vector<path_handle_t> &,
                               const xp::XP &)> get_sum_path_lengths
                = [&](const std::vector<path_handle_t> &path_sgd_use_paths, const xp::XP &path_index) {
                    uint64_t sum_path_length = 0;
                    for (auto &path : path_sgd_use_paths) {
                        sum_path_length += path_index.get_path_length(path);
                    }
                    return sum_path_length;
                };
        std::function<uint64_t(const std::vector<path_handle_t> &,
                               const xp::XP &)> get_max_path_length
                = [&](const std::vector<path_handle_t> &path_sgd_use_paths, const xp::XP &path_index) {
                    uint64_t max_path_length = 0;
                    for (auto &path : path_sgd_use_paths) {
                        max_path_length = std::max(max_path_length, path_index.get_path_length(path));
                    }
                    return max_path_length;
                };
        // default parameters
        std::string path_sgd_seed;
        if (p_sgd_seed) {
            if (num_threads > 1) {
                std::cerr
                        << "[odgi layout2] Error: Please only specify a seed for the path guided 2D linear SGD when using 1 thread."
                        << std::endl;
                return 1;
            }
            path_sgd_seed = args::get(p_sgd_seed);
        } else {
            path_sgd_seed = "pangenomic!";
        }
        if (p_sgd_min_term_updates_paths && p_sgd_min_term_updates_num_nodes) {
            std::cerr
                    << "[odgi layout2] Error: There can only be one argument provided for the minimum number of term updates in the path guided 1D SGD."
                       "Please either use -G=[N], path-sgd-min-term-updates-paths=[N] or -U=[N], path-sgd-min-term-updates-nodes=[N]."
                    << std::endl;
            return 1;
        }
        uint64_t path_sgd_iter_max = args::get(p_sgd_iter_max) ? args::get(p_sgd_iter_max) : 30;
        uint64_t path_sgd_iter_max_learning_rate = args::get(p_sgd_iter_with_max_learning_rate) ? args::get(
                p_sgd_iter_with_max_learning_rate) : 0;
        double path_sgd_zipf_theta = args::get(p_sgd_zipf_theta) ? args::get(p_sgd_zipf_theta) : 0.99;
        double path_sgd_eps = args::get(p_sgd_eps) ? args::get(p_sgd_eps) : 0.01;
        double path_sgd_delta = args::get(p_sgd_delta) ? args::get(p_sgd_delta) : 0;
        // will be filled, if the user decides to write a snapshot of the graph after each sorting iterationn
        std::vector<std::vector<handle_t>> snapshots;
        const bool snapshot = p_sgd_snapshot;
        const bool path_sgd_deterministic = p_sgd_deterministic;

        // default parameters that need a path index to be present
        uint64_t path_sgd_min_term_updates;
        uint64_t path_sgd_zipf_space;
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
                    std::cerr << "[odgi layout2] Error: Path '" << buf
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
        uint64_t max_path_length = get_max_path_length(path_sgd_use_paths, path_index);
        path_sgd_zipf_space = args::get(p_sgd_zipf_space) ? args::get(p_sgd_zipf_space) : max_path_length;
        double path_sgd_max_eta = args::get(p_sgd_eta_max) ? args::get(p_sgd_eta_max) : max_path_length * max_path_length;


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
        // std::uniform_int_distribution<uint64_t> dist(1, graph.get_node_count());
        std::normal_distribution<double> dist(0,1); //graph.get_node_count());

        uint64_t len = 0;
        nid_t last_node_id = graph.min_node_id();
        graph.for_each_handle([&](const handle_t &h) {
            nid_t node_id = graph.get_id(h);
            if (node_id - last_node_id > 1) {
                std::cerr << "[odgi layout2] error: The graph is not optimized. Please run 'odgi sort' using -O, --optimize" << std::endl;
                exit(1);
            }
            last_node_id = node_id;

            uint64_t pos = 2 * number_bool_packing::unpack_number(h);
            //double x_offset = len;
            //double y_offset = dist(rng);
            graph_X[pos].store(len); //dist(rng));
            graph_Y[pos].store(dist(rng));
            len += graph.get_length(h);
            graph_X[pos + 1].store(len); //dist(rng));
            graph_Y[pos + 1].store(dist(rng));

            //std::cerr << pos << ": " << graph_X[pos] << "," << graph_Y[pos] << " ------ " << graph_X[pos + 1] << "," << graph_Y[pos + 1] << std::endl;
        });

        //double max_x = 0;

        std::vector<std::vector<double>> snapshotsX; // TODO to remove
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
                num_threads,
                true,
                snapshot,
                snapshotsX,
                graph_X,
                graph_Y
            );

        /*for (auto &component_order : weak_component_order) {
            auto &weak_component = weak_components[component_order.second];

            std::vector<handlegraph::nid_t> component_ids;
            for (auto &id : weak_component) {
                component_ids.push_back(id);
            }

            std::sort(component_ids.begin(), component_ids.end());
            const handle_t &first_handle = graph.get_handle(component_ids.front());
            uint64_t pos_offset = number_bool_packing::unpack_number(first_handle);

            uint64_t n = weak_component.size();

            // Generate randomly X and Y
            std::vector<double> component_X(2 * n);
            std::vector<double> component_Y(2 * n);

            uint64_t len = 0;
            for (auto node_id : weak_component) {
                const handle_t &h = graph.get_handle(node_id);

                uint64_t pos = number_bool_packing::unpack_number(h) - pos_offset;

                graph_X[pos] = len;     // node+
                graph_Y[pos] = dist(rng);
                len += graph.get_length(h);
                graph_X[pos + 1] = len; // node-
                graph_Y[pos + 1] = dist(rng); // node-
            }

            ///std::vector<double> component_layout = algorithms::path_2D_sgd_order(graph,
            //                                                  path_index,
            //                                                  path_sgd_use_paths,
            //                                                  path_sgd_iter_max,
            //                                                  path_sgd_iter_max_learning_rate,
            //                                                  path_sgd_min_term_updates,
            //                                                  path_sgd_delta,
            //                                                  path_sgd_eps,
            //                                                  path_sgd_max_eta,
            //                                                  path_sgd_zipf_theta,
            //                                                  path_sgd_zipf_space,
            //                                                  num_threads,
            //                                                  progress,
            //                                                  path_sgd_seed,
            //                                                  snapshot,
            //                                                  snapshots,
            //                                                  component_X,
            //                                                  component_Y);


            // ToDo Formula: 2(id-1) + alfa (if reverse)
            for (uint64_t i = 0; i < 2 * n; i += 2) {
                //std::cerr << "i = " << i << std::endl;
                uint64_t j = component_ids[i / 2] - 1;

                graph_X[j * 2] = component_X[i] + max_x;
                graph_X[j * 2 + 1] = component_X[i] + max_x;

                graph_Y[j * 2] = component_Y[i];
                graph_Y[j * 2 + 1] = component_Y[i];

                //std::cerr << "layout " << j << " " << layout[j*2] << " " << layout[j*2+1] << std::endl;
            }
            // set new max_x
            for (uint64_t i = 0; i < 2 * n; i += 2) {
                max_x = std::max(component_X[i], max_x);
            }
            max_x += x_padding;
        }*/



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

        std::string outfile = args::get(svg_out_file);
        if (outfile.size()) {
            if (outfile == "-") {
                draw_svg(std::cout, X_final, Y_final, graph, svg_scale);
            } else {
                ofstream f(outfile.c_str());
                //draw_svg(f, X_final, Y_final, graph, svg_scale);
                to_tsv(f, X_final, Y_final, graph, svg_scale);
                f.close();
            }
        }
        return 0;
    }

    static Subcommand odgi_layout2("layout2", "use SGD to make 2D layouts of the graph",
                                   PIPELINE, 3, main_layout2);


}
