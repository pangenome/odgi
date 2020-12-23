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
    args::ValueFlag<std::string> mean_dists_out_file(parser, "FILE", "write the mean distances of adjacent nodes for each node to this file", {'M', "mean-dists-adj-nodes"});
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

    if (mean_dists_out_file) {
        std::string distance_out = args::get(mean_dists_out_file);

        /// distances of consecutive nodes, will work for both 1D and 2D!
        vector<double> consecutive_node_dists;
        uint64_t handle_pos = 0;

        /// average distances to all adjacent nodes for each node
        vector<double> mean_dists_adj_nodes(graph.get_node_count(), 0.0);

        graph.for_each_handle(
                [&](const handle_t &h) {
                    vector<double> adj_node_dists;
                    // we didn't hit the last node, yet
                    if (handle_pos < graph.max_node_id() - 1) {
                        double dx = layout.get_x(handle_pos) - layout.get_x(handle_pos + 1);
                        double dy = layout.get_y(handle_pos) - layout.get_y(handle_pos + 1);
                        double handle_dist = sqrt(dx * dx + dy * dy);
                        consecutive_node_dists.push_back(handle_dist);
                    }
                    // collect distances to the left
                    graph.follow_edges(h, true, [&](const handle_t &o) {
                        nid_t o_id = graph.get_id(o);
                        double dx = layout.get_x(handle_pos) - layout.get_x(o_id - 1);
                        double dy = layout.get_y(handle_pos) - layout.get_y(o_id - 1);
                        double handle_dist = sqrt(dx * dx + dy * dy);
                        adj_node_dists.push_back(handle_dist);
                    });
                    // collect distances to the right
                    graph.follow_edges(h, false, [&](const handle_t &o) {
                        nid_t o_id = graph.get_id(o);
                        double dx = layout.get_x(handle_pos) - layout.get_x(o_id - 1);
                        double dy = layout.get_y(handle_pos) - layout.get_y(o_id - 1);
                        double handle_dist = sqrt(dx * dx + dy * dy);
                        adj_node_dists.push_back(handle_dist);
                    });
                    // calculate mean distance and add it to the overall distances
                    double sum = std::accumulate(adj_node_dists.begin(), adj_node_dists.end(), 0.0);
                    double mean = sum / adj_node_dists.size();
                    mean_dists_adj_nodes[handle_pos] = mean;
                    handle_pos++;
                });

        /// calculate average and median distances of consecutive nodes and print
        double sum = std::accumulate(consecutive_node_dists.begin(), consecutive_node_dists.end(), 0.0);
        double mean = sum / consecutive_node_dists.size();
        std::cerr << "[odgi cut]: average distance of consecutive nodes: " << mean << std::endl;
        const auto median_it = consecutive_node_dists.begin() + consecutive_node_dists.size() / 2;
        std::nth_element(consecutive_node_dists.begin(), median_it, consecutive_node_dists.end());
        auto median = *median_it;
        std::cerr << "[odgi cut]: median distance of consecutive nodes: " << median << std::endl;

        ofstream f(distance_out.c_str());
        uint64_t index = 0;
        f << "idx" << "\t" << "mean_dist_adj_nodes" << std::endl;
        for (auto &dist : mean_dists_adj_nodes) {
            f << index << "\t" << dist << std::endl;
            index++;
        }
        f.close();
    }

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
