#include "subcommand.hpp"
#include "graph.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "algorithms/hash.hpp"
#include "lodepng.h"
#include <iostream>

namespace odgi {

namespace png {

/*
3 ways to encode a PNG from RGBA pixel data to a file (and 2 in-memory ways).
NOTE: this samples overwrite the file or test.png without warning!
*/

//g++ lodepng.cpp examples/example_encode.cpp -I./ -ansi -pedantic -Wall -Wextra -O3

//Example 1
//Encode from raw pixels to disk with a single function call
//The image argument has width * height RGBA pixels or width * height * 4 bytes
void encodeOneStep(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
  //Encode the image
  unsigned error = lodepng::encode(filename, image, width, height);

  //if there's an error, display it
  if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

//Example 2
//Encode from raw pixels to an in-memory PNG file first, then write it to disk
//The image argument has width * height RGBA pixels or width * height * 4 bytes
void encodeTwoSteps(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
  std::vector<unsigned char> png;

  unsigned error = lodepng::encode(png, image, width, height);
  if(!error) lodepng::save_file(png, filename);

  //if there's an error, display it
  if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

//Example 3
//Save a PNG file to disk using a State, normally needed for more advanced usage.
//The image argument has width * height RGBA pixels or width * height * 4 bytes
void encodeWithState(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
  std::vector<unsigned char> png;
  lodepng::State state; //optionally customize this one

  unsigned error = lodepng::encode(png, image, width, height, state);
  if(!error) lodepng::save_file(png, filename);

  //if there's an error, display it
  if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

}


using namespace odgi::subcommand;

int main_viz(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi viz";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("variation graph visualizations");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the index from this file", {'i', "idx"});
    args::ValueFlag<std::string> png_out_file(parser, "FILE", "write the output (png) to this file", {'o', "out"});
    args::ValueFlag<uint64_t> image_width(parser, "N", "width in pixels of output image", {'x', "width"});
    args::ValueFlag<uint64_t> image_height(parser, "N", "height in pixels of output image", {'y', "height"});
    args::ValueFlag<uint64_t> path_height(parser, "N", "path display height", {'P', "path-height"});
    args::ValueFlag<uint64_t> threads(parser, "N", "number of threads to use", {'t', "threads"});

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

    if (args::get(threads)) {
        omp_set_num_threads(args::get(threads));
    } else {
        omp_set_num_threads(1);
    }

    graph_t graph;
    assert(argc > 0);
    std::string infile = args::get(dg_in_file);
    if (infile.size()) {
        ifstream f(infile.c_str());
        graph.load(f);
        f.close();
    }

    //NOTE: this sample will overwrite the file or test.png without warning!
    //const char* filename = argc > 1 ? argv[1] : "test.png";
    if (!args::get(png_out_file).size()) {
        std::cerr << "[odgi viz] error: output image required" << std::endl;
        return 1;
    }
    const char* filename = args::get(png_out_file).c_str();

    std::vector<uint64_t> position_map(graph.node_size()+1);
    std::vector<std::pair<uint64_t, uint64_t>> contacts;
    uint64_t len = 0;
    graph.for_each_handle([&](const handle_t& h) {
            position_map[number_bool_packing::unpack_number(h)] = len;
            uint64_t hl = graph.get_length(h);
            len += hl;
        });
    position_map[position_map.size()-1] = len;

    uint64_t path_count = graph.get_path_count();
    uint64_t pix_per_path = (args::get(path_height) ? args::get(path_height) : 10);
    uint64_t path_space = path_count * pix_per_path;
    // the math here only works if the image size is greater than or equal to the graph length
    uint64_t width = std::min(len, (args::get(image_width) ? args::get(image_width) : 1000));
    uint64_t height = std::min(len, (args::get(image_height) ? args::get(image_height) : 1000));
    std::vector<uint8_t> image;
    image.resize(width * (height + path_space) * 4, 255);
    float scale_x = (float)width/(float)len;
    float scale_y = (float)height/(float)len;

    auto add_point = [&](const uint64_t& _x, const uint64_t& _y,
                         const uint8_t& _r, const uint8_t& _g, const uint8_t& _b) {
        uint64_t x = std::min((uint64_t)std::round(_x * scale_x), width-1);
        uint64_t y = std::min((uint64_t)std::round(_y * scale_y), height-1) + path_space;
        image[4 * width * y + 4 * x + 0] = _r;
        image[4 * width * y + 4 * x + 1] = _g;
        image[4 * width * y + 4 * x + 2] = _b;
        image[4 * width * y + 4 * x + 3] = 255;
    };

    auto add_edge = [&](const handle_t& h, const handle_t& o) {
        auto& _a = position_map[number_bool_packing::unpack_number(h) + !number_bool_packing::unpack_bit(h)];
        auto& _b = position_map[number_bool_packing::unpack_number(o) + !number_bool_packing::unpack_bit(o)];
        uint64_t a = std::min(_a, _b);
        uint64_t b = std::max(_a, _b);
        uint64_t dist = b - a;
        uint64_t i = 0;
        for ( ; i < dist; i+=1/scale_y) {
            add_point(a, i, 0, 0, 0);
        }
        while (a < b) {
            add_point(a, i, 0, 0, 0);
            a += 1/scale_x;
        }
        for (uint64_t j = 0; j < dist; j+=1/scale_y) {
            add_point(b, j, 0, 0, 0);
        }
    };

    graph.for_each_handle([&](const handle_t& h) {
            uint64_t p = position_map[number_bool_packing::unpack_number(h)];
            uint64_t hl = graph.get_length(h);
            // make contects for the bases in the node
            //for (uint64_t i = 0; i < hl; ++i) {
            for (uint64_t i = 0; i < hl; i+=1/scale_x) {
                add_point(p+i, 0, 0, 0, 0);
            }
        });

    graph.for_each_handle([&](const handle_t& h) {
            // add contacts for the edges
            graph.follow_edges(h, false, [&](const handle_t& o) {
                    add_edge(h, o);
                });
            graph.follow_edges(h, true, [&](const handle_t& o) {
                    add_edge(o, h);
                });
        });

    auto add_path_step = [&](const uint64_t& _x, const uint64_t& _y,
                             const uint8_t& _r, const uint8_t& _g, const uint8_t& _b) {
        uint64_t x = std::min((uint64_t)std::round(_x * scale_x), width-1);
        uint64_t t = _y*pix_per_path;
        uint64_t s = t+pix_per_path;
        for (uint64_t y = t; y < s; ++y) {
            image[4 * width * y + 4 * x + 0] = _r;
            image[4 * width * y + 4 * x + 1] = _g;
            image[4 * width * y + 4 * x + 2] = _b;
            image[4 * width * y + 4 * x + 3] = 255;
        }
    };
    
    graph.for_each_path_handle([&](const path_handle_t& path) {
            std::hash<std::string> hash_fn;
            uint64_t path_name_hash_r = wang_hash_64(hash_fn(graph.get_path_name(path)+"red"));
            uint64_t path_name_hash_g = wang_hash_64(hash_fn(graph.get_path_name(path)+"green"));
            uint64_t path_name_hash_b = wang_hash_64(hash_fn(graph.get_path_name(path)+"blue"));
            float path_r_f = (float)path_name_hash_r/(float)(std::numeric_limits<uint64_t>::max());
            float path_g_f = (float)path_name_hash_g/(float)(std::numeric_limits<uint64_t>::max());
            float path_b_f = (float)path_name_hash_b/(float)(std::numeric_limits<uint64_t>::max());
            float sum = path_r_f + path_g_f + path_b_f;
            path_r_f /= sum;
            path_g_f /= sum;
            path_b_f /= sum;
            // brighten the color
            float f = std::min(1.8, 1.0/std::max(std::max(path_r_f, path_g_f), path_b_f));
            uint8_t path_r = (uint8_t)std::round(255*std::min(path_r_f*f, (float)1.0));
            uint8_t path_g = (uint8_t)std::round(255*std::min(path_g_f*f, (float)1.0));
            uint8_t path_b = (uint8_t)std::round(255*std::min(path_b_f*f, (float)1.0));
            std::cerr << "path " << as_integer(path) << " " << graph.get_path_name(path) << " " << path_r_f << " " << path_g_f << " " << path_b_f
                      << " " << (int)path_r << " " << (int)path_g << " " << (int)path_b << std::endl;
            /// Loop over all the occurrences along a path, from first through last
            uint64_t path_rank = as_integer(path);
            uint64_t step = 0;
            graph.for_each_occurrence_in_path(path, [&](const occurrence_handle_t& occ) {
                    handle_t h = graph.get_occurrence(occ);
                    uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                    // todo...
                    uint64_t hl = graph.get_length(h);
                    // make contects for the bases in the node
                    for (uint64_t i = 0; i < hl; i+=1/scale_x) {
                        add_path_step(p+i, path_rank, path_r, path_g, path_b);
                    }
                 });
         });
    
    png::encodeOneStep(filename, image, width, height + path_space);

    return 0;
}

static Subcommand odgi_viz("viz", "visualize the graph",
                           PIPELINE, 3, main_viz);


}
