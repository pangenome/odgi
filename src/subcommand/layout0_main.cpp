#include "subcommand.hpp"
#include <iostream>
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/sgd_layout.hpp"

namespace odgi {

using namespace odgi::subcommand;

void draw_svg(ostream& out, const std::vector<double>& X, const HandleGraph& graph, double scale) {
    double border = 10.0;
    double min_x = 0;
    double min_y = 0;
    double max_x = 0;
    double max_y = 0;
    // determine boundaries
    uint64_t n = graph.get_node_count();
    for (uint64_t i = 0; i < n; ++i) {
        double x = X[i*2]*scale;
        double y = X[i*2+1]*scale;
        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
    }
    double width = max_x - min_x;
    double height = max_y - min_y;
    
    out << "<svg width=\"" << width + border << "\" height=\"" << height + border << "\" "
        << "viewBox=\"" << min_x - border/2<< " " << min_y - border/2 << " " << width + border << " " << height + border << "\" xmlns=\"http://www.w3.org/2000/svg\">"
        << "<style type=\"text/css\">"
        << "line{stroke:black;stroke-width:1.0;stroke-opacity:1.0;stroke-linecap:round;}"
        //<< "circle{{r:" << 1.0 << ";fill:black;fill-opacity:" << 1.0 << "}}"
        << "</style>"
        << std::endl;

    // iterate through graph edges
    graph.for_each_edge([&](const edge_t& e) {
            uint64_t a = graph.get_id(e.first)-1;
            uint64_t b = graph.get_id(e.second)-1;
            out << "<line x1=\"" << X[a*2]*scale << "\" x2=\"" << X[b*2]*scale
                << "\" y1=\"" << X[a*2+1]*scale << "\" y2=\"" << X[b*2+1]*scale << "\"/>"
                << std::endl;
        });
    /* // to draw nodes
    for (uint64_t i = 0; i < n; ++i) {
        std::cout << "<circle cx=\"" << X[i*2]*scale << "\" cy=\"" << X[i*2+1]*scale << "\" r=\"1.0\"/>" << std::endl;
    }
    */
    out << "</svg>" << std::endl;
}

int main_layout0(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi layout0";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("draw 2D layouts of the graph using stochastic gradient descent (the graph must be sorted and id-compacted)");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> svg_out_file(parser, "FILE", "write the SVG rendering to this file", {'o', "out"});
    args::ValueFlag<uint64_t> iter_max(parser, "N", "maximum number of iterations to run the layout (default 30)", {'m', "iter-max"});
    args::ValueFlag<uint64_t> num_pivots(parser, "N", "number of pivots for sparse layout (default 0, non-sparse layout)", {'p', "n-pivots"});
    args::ValueFlag<double> eps_rate(parser, "N", "learning rate for SGD layout (default 0.01)", {'e', "eps"});
    args::ValueFlag<double> x_pad(parser, "N", "padding between connected component layouts (default 10.0)", {'x', "x-padding"});
    args::ValueFlag<double> render_scale(parser, "N", "SVG scaling (default 5.0)", {'R', "render-scale"});
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
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    if (!dg_in_file) {
        std::cerr << "[odgi::layout0] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (!svg_out_file) {
        std::cerr << "[odgi::layout0] error: please specify an output file to where to store the layout via -o=[FILE], --out=[FILE]." << std::endl;
        return 1;
    }

    uint64_t t_max = !args::get(iter_max) ? 30 : args::get(iter_max);
    uint64_t n_pivots = !args::get(num_pivots) ? 0 : args::get(num_pivots);
    double eps = !args::get(eps_rate) ? 0.01 : args::get(eps_rate);
    double x_padding = !args::get(x_pad) ? 10.0 : args::get(x_pad);
    double svg_scale = !args::get(render_scale) ? 5.0 : args::get(render_scale);
    
    graph_t graph;
    assert(argc > 0);
    std::string infile = args::get(dg_in_file);
    if (!infile.empty()) {
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
            ifstream f(infile.c_str());
            graph.deserialize(f);
            f.close();
        }
    }

    std::vector<double> layout = algorithms::sgd_layout(graph, n_pivots, t_max, eps, x_padding);

    std::string outfile = args::get(svg_out_file);
    if (!outfile.empty()) {
        if (outfile == "-") {
            draw_svg(std::cout, layout, graph, svg_scale);
        } else {
            ofstream f(outfile.c_str());
            draw_svg(f, layout, graph, svg_scale);
            f.close();
        }
    }
    return 0;
}

static Subcommand odgi_layout0("layout0", "use SGD to make 2D layouts of the graph",
                               DEVELOPMENT, 3, main_layout0);


}
