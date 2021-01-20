#include "subcommand.hpp"
#include <iostream>
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/layout.hpp"
#include "algorithms/bed_records.hpp"
#include <numeric>

namespace odgi {

using namespace odgi::subcommand;

int main_tension(int argc, char **argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc - 1; ++i) {
        argv[i] = argv[i + 1];
    }
    std::string prog_name = "odgi tension";
    argv[0] = (char *) prog_name.c_str();
    --argc;

    args::ArgumentParser parser(
        "evaluate the tension of a graph helping to locate structural variants and abnormalities");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> layout_in_file(parser, "FILE", "read the layout coordinates from this .lay format file produced by odgi sort or odgi layout", {'c', "coords-in"});
    args::ValueFlag<double> window_size(parser, "N", "window size in kb in which each tension is calculated", {'w', "window-size"});
    args::ValueFlag<std::string> tsv_out_file(parser, "FILE", "write the TSV layout to this file", {'T', "tsv"});
    // args::ValueFlag<std::string> bed_out_file(parser, "FILE", "write the BED intervals to this file", {'B', "bed"}); we write to stdout
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

    uint64_t thread_count = 1;
    if (args::get(nthreads)) {
        omp_set_num_threads(args::get(nthreads));
        thread_count = args::get(nthreads);
    }

    if (!dg_in_file) {
        std::cerr
                << "[odgi tension] error: Please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
                << std::endl;
        return 1;
    }

    if (!window_size) {
        std::cerr
                << "[odgi tension] error: Please specify a window size for the detection of the tension level of the graph via -w=[N], --window-size=[N]."
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

    double window_size_ = args::get(window_size) * 1000;

    vector<path_handle_t> p_handles;
    graph.for_each_path_handle([&] (const path_handle_t &p) {
       p_handles.push_back(p);
    });

    // algorithms::bed_records_class bed;

    // bed.open_writer();
// #pragma omp parallel for schedule(static, 1) num_threads(thread_count)
    for (uint64_t idx = 0; idx < p_handles.size(); idx++) {
        path_handle_t p = p_handles[idx];
        std::string path_name = graph.get_path_name(p);
        uint64_t cur_window_start = 1;
        uint64_t cur_window_end = 0;
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
                uint64_t nuc_dist = graph.get_length(h);
                path_nuc_dist += nuc_dist;
                cur_window_end += nuc_dist;
                if ((cur_window_end - cur_window_start + 1) >= window_size_) {
                    // BED files are 0-based http://genome.ucsc.edu/FAQ/FAQformat#format1
                    std::cout << path_name << "\t" // chrom
                              << (cur_window_start - 1) << "\t" // chromStart
                              << (cur_window_end) << "\t" // chromEnd
                              // << path_name << ":" << path_layout_dist << "-" << path_nuc_dist << "\t" // name
                              // << 1000 << "\t" // score
                              //<< "+" << "\t" // strand
                              << path_layout_dist << "\t"
                              << path_nuc_dist
                              << std::endl;
                    cur_window_start = cur_window_end + 1;
                    cur_window_end = cur_window_start - 1;
                    path_layout_dist = 0;
                    path_nuc_dist = 0;
                }
            } else {
                // we only take a look at the current node
                /// f
                if (!graph.get_is_reverse(h)) {
                    path_layout_dist += algorithms::layout::coord_dist(h_coords_start, h_coords_end);
                } else {
                    /// r
                    path_layout_dist += algorithms::layout::coord_dist(h_coords_end, h_coords_start);
                }
                uint64_t nuc_dist = graph.get_length(h);
                path_nuc_dist += nuc_dist;
                cur_window_end += nuc_dist;
                // BED files are 0-based http://genome.ucsc.edu/FAQ/FAQformat#format1
                std::cout << path_name << "\t" // chrom
                          << (cur_window_start - 1) << "\t" // chromStart
                          << (cur_window_end) << "\t" // chromEnd
                          // << path_name << ":" << path_layout_dist << "-" << path_nuc_dist << "\t" // name
                          // << 1000 << "\t" // score
                          //<< "+" << "\t" // strand
                          << path_layout_dist << "\t"
                          << path_nuc_dist
                          << std::endl;
            }
        });
    }
    // bed.close_writer();

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

static Subcommand odgi_tension("tension", "evaluate the tension of a graph helping to locate structural variants and abnormalities",
                            PIPELINE, 3, main_tension);


}
