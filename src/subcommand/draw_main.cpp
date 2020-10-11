#include "subcommand.hpp"
#include <iostream>
#include "odgi.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "algorithms/xp.hpp"
#include "algorithms/draw.hpp"
#include "algorithms/layout.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_draw(int argc, char **argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc - 1; ++i) {
        argv[i] = argv[i + 1];
    }
    std::string prog_name = "odgi draw";
    argv[0] = (char *) prog_name.c_str();
    --argc;

    args::ArgumentParser parser(
        "draw previously-determined 2D layouts of the graph with diverse annotations");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> layout_in_file(parser, "FILE", "read the layout coordinates from this .lay format file produced by odgi layout", {'c', "coords-in"});
    //args::Flag in_is_tsv(parser, "is-tsv", "if the input is .tsv format (three column: id, X, Y) rather the default .lay binary format", {'I', "input-is-tsv"});
    args::ValueFlag<std::string> tsv_out_file(parser, "FILE", "write the TSV layout plus displayed annotations to this file", {'T', "tsv"});
    args::ValueFlag<std::string> svg_out_file(parser, "FILE", "write an SVG rendering to this file", {'s', "svg"});
    args::ValueFlag<std::string> png_out_file(parser, "FILE", "write a rasterized PNG rendering to this file", {'p', "png"});
    args::ValueFlag<uint64_t> png_height(parser, "FILE", "height of PNG rendering (default: 1000)", {'H', "png-height"});
    args::ValueFlag<uint64_t> png_border(parser, "FILE", "size of PNG border in bp (default: 10)", {'E', "png-border"});
    args::Flag color_paths(parser, "color-paths", "color paths (in PNG output)", {'C', "color-paths"});
    args::ValueFlag<double> render_scale(parser, "N", "image scaling (default 1.0)", {'R', "scale"});
    args::ValueFlag<double> render_border(parser, "N", "image border (in approximate bp) (default 100.0)", {'B', "border"});
    args::ValueFlag<double> png_line_width(parser, "N", "line width (in approximate bp) (default 0.0)", {'w', "line-width"});
    args::ValueFlag<double> png_line_overlay(parser, "N", "line width (in approximate bp) (default 10.0)", {'O', "line-overlay"});
    args::ValueFlag<double> png_path_line_spacing(parser, "N", "spacing between path lines in png layout (in approximate bp) (default 1.0)", {'S', "path-line-spacing"});
    args::ValueFlag<std::string> xp_in_file(parser, "FILE", "load the path index from this file", {'X', "path-index"});
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
            << "[odgi draw] error: Please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
            << std::endl;
        return 1;
    }

    if (!tsv_out_file && !svg_out_file && !png_out_file) {
        std::cerr
            << "[odgi draw] error: Please specify an output file to where to store the layout via -p/--png=[FILE], -s/--svg=[FILE], -T/--tsv=[FILE]"
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

    const double svg_scale = !render_scale ? 1.0 : args::get(render_scale);
    const double border_bp = !render_border ? 100.0 : args::get(render_border);

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

    if (svg_out_file) {
        auto& outfile = args::get(svg_out_file);
        ofstream f(outfile.c_str());
        // todo could be done with callbacks
        std::vector<double> X = layout.get_X();
        std::vector<double> Y = layout.get_Y();
        algorithms::draw_svg(f, X, Y, graph, svg_scale, border_bp);
        f.close();    
    }

    if (png_out_file) {
        auto& outfile = args::get(png_out_file);
        uint64_t _png_height = png_height ? args::get(png_height) : 1000;
        double _png_line_width = png_line_width ? args::get(png_line_width) : 0;
        bool _color_paths = args::get(color_paths);
        double _png_path_line_spacing = png_path_line_spacing ? args::get(png_path_line_spacing) : 1.0;
        double _png_line_overlay = png_line_overlay ? args::get(png_line_overlay) : 10.0;
        // todo could be done with callbacks
        std::vector<double> X = layout.get_X();
        std::vector<double> Y = layout.get_Y();
        algorithms::draw_png(outfile, X, Y, graph, 1.0, border_bp, 0, _png_height, _png_line_width, _png_line_overlay, _png_path_line_spacing, _color_paths);
    }
    
    return 0;
}

static Subcommand odgi_draw("draw", "draw previously-determined 2D layouts of the graph with diverse annotations",
                            PIPELINE, 3, main_draw);


}
