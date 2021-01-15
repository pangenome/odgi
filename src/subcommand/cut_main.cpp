#include "subcommand.hpp"
#include <iostream>
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/layout.hpp"
#include <numeric>

namespace odgi {

using namespace odgi::subcommand;

int main_cut(int argc, char **argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc - 1; ++i) {
        argv[i] = argv[i + 1];
    }
    std::string prog_name = "odgi cut";
    argv[0] = (char *) prog_name.c_str();
    --argc;

    args::ArgumentParser parser(
        "cut out structural variants by given distance threshold");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> layout_in_file(parser, "FILE", "read the layout coordinates from this .lay format file produced by odgi sort or odgi layout", {'c', "coords-in"});
    args::ValueFlag<std::string> tsv_out_file(parser, "FILE", "write the TSV layout to this file", {'T', "tsv"});
    args::Flag progress(parser, "progress", "display progress of the sort", {'P', "progress"});
    args::ValueFlag<uint64_t> nthreads(parser, "N", "number of threads to use for parallel phases", {'t', "threads"});
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

    if (args::get(nthreads)) {
        omp_set_num_threads(args::get(nthreads));
    }

    if (!dg_in_file) {
        std::cerr
            << "[odgi cut] error: Please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
            << std::endl;
        return 1;
    }

    if (!tsv_out_file) {
        std::cerr
            << "[odgi cut] error: Please specify an output file to where to store the layout via -p/--png=[FILE], -s/--svg=[FILE], -T/--tsv=[FILE]"
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

    algorithms::layout::Layout layout;
    if (layout_in_file) {
        auto& infile = args::get(layout_in_file);
        if (infile.size()) {
            if (infile == "-") {
                layout.load(std::cin);
            } else {
                ifstream f(infile.c_str());
                layout.load(f);
                f.close();
            }
        }
    }

    std::cout << "path_name\tlayout_distance\tnucleotide_distance" << std::endl;

    graph.for_each_path_handle([&](const path_handle_t &p) {
        std::string path_name = graph.get_path_name(p);
        double path_layout_dist = 0;
        uint64_t path_nuc_dist = 0;
        graph.for_each_step_in_path(p, [&](const step_handle_t &s) {
            handle_t h = graph.get_handle_of_step(s);
            algorithms::xy_d_t h_coords_start;
            algorithms::xy_d_t h_coords_end;
            if (graph.get_is_reverse(h)) {
                h_coords_start = layout.coords(graph.flip(h));
                h_coords_end = layout.coords(h);
            } else {
                h_coords_start = layout.coords(h);
                h_coords_end = layout.coords(graph.flip(h));
            }
            // did we hit the last step?
            if (graph.has_next_step(s)) {
                step_handle_t next_s = graph.get_next_step(s);
                handle_t next_h = graph.get_handle_of_step(next_s);
                algorithms::xy_d_t next_h_coords_start;
                algorithms::xy_d_t next_h_coords_end;
                if (graph.get_is_reverse(next_h)) {
                    next_h_coords_start = layout.coords(graph.flip(next_h));
                    next_h_coords_end = layout.coords(next_h);
                } else {
                    next_h_coords_start = layout.coords(next_h);
                    next_h_coords_end = layout.coords(graph.flip(next_h));
                }
                double within_node_dist = 0;
                double from_node_to_node_dist = 0;
                if (!graph.get_is_reverse(h)) {
                    /// f + f
                    if (!graph.get_is_reverse(next_h)) {
                        within_node_dist = algorithms::layout::coord_dist(h_coords_start, h_coords_end);
                        from_node_to_node_dist = algorithms::layout::coord_dist(h_coords_end, next_h_coords_start);
                    } else {
                        /// f + r
                        within_node_dist = algorithms::layout::coord_dist(h_coords_start, h_coords_end);
                        from_node_to_node_dist = algorithms::layout::coord_dist(h_coords_end, next_h_coords_end);
                    }
                } else {
                    /// r + r
                    if (graph.get_is_reverse(next_h)) {
                        within_node_dist = algorithms::layout::coord_dist(h_coords_end, h_coords_start);
                        from_node_to_node_dist = algorithms::layout::coord_dist(h_coords_start, next_h_coords_start);
                    } else {
                        /// r + f
                        within_node_dist = algorithms::layout::coord_dist(h_coords_end, h_coords_start);
                        from_node_to_node_dist = algorithms::layout::coord_dist(h_coords_start, next_h_coords_end);
                    }
                }
                path_layout_dist += within_node_dist;
                path_layout_dist += from_node_to_node_dist;
                path_nuc_dist += graph.get_length(h);
            } else {
                // we only take a look at the current node
                /// f
                if (!graph.get_is_reverse(h)) {
                    path_layout_dist += algorithms::layout::coord_dist(h_coords_start, h_coords_end);
                } else {
                    /// r
                    path_layout_dist += algorithms::layout::coord_dist(h_coords_end, h_coords_start);
                }
                path_nuc_dist += graph.get_length(h);
            }
        });
        std::cout << path_name << "\t" << path_layout_dist << "\t" << path_nuc_dist << std::endl;
    });

    if (tsv_out_file) {
        auto& outfile = args::get(tsv_out_file);
        if (outfile.size()) {
            if (outfile == "-") {
                layout.to_tsv(std::cout);
            } else {
                ofstream f(outfile.c_str());
                layout.to_tsv(f);
                f.close();
            }
        }
    }

    return 0;
}

static Subcommand odgi_draw("cut", "cut out structural variants by given distance threshold",
                            PIPELINE, 3, main_cut);


}
