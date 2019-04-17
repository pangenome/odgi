#include "subcommand.hpp"
#include "graph.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "algorithms/hash.hpp"
#include "algorithms/id_ordered_paths.hpp"
#include "lodepng.h"
#include <limits>
//#include "picosha2.h"
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
    args::ValueFlag<uint64_t> path_x_pad(parser, "N", "path x padding", {'X', "path-x-padding"});
    args::Flag path_per_row(parser, "bool", "display a single path per row rather than packing them", {'R', "path-per-row"});
    args::ValueFlag<float> link_path_pieces(parser, "FLOAT", "show thin links of this relative width to connect path pieces", {'L', "link-path-pieces"});
    args::ValueFlag<std::string> alignment_prefix(parser, "STRING", "apply alignment-related visual motifs to paths with this name prefix", {'A', "alignment-prefix"});
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
    uint64_t pix_per_link = std::max((uint64_t)1, (uint64_t)std::round(args::get(link_path_pieces)*pix_per_path));
    uint64_t link_pix_y = pix_per_path/2 - pix_per_link/2;
    uint64_t path_space = path_count * pix_per_path;
    uint64_t path_padding = args::get(path_x_pad);
    uint64_t bottom_padding = 5;
    // the math here only works if the image size is greater than or equal to the graph length
    uint64_t width = std::min(len, (args::get(image_width) ? args::get(image_width) : 1000));
    uint64_t height = std::min(len, (args::get(image_height) ? args::get(image_height) : 1000) + bottom_padding);
    std::vector<uint8_t> image;
    image.resize(width * (height + path_space) * 4, 255);
    float scale_x = (float)width/(float)len;
    float scale_y = (float)height/(float)len;

    bool aln_mode = !args::get(alignment_prefix).empty();
    std::string aln_prefix;
    if (aln_mode) {
        aln_prefix = args::get(alignment_prefix);
    }

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

    auto add_path_link = [&](const uint64_t& _x, const uint64_t& _y,
                             const uint8_t& _r, const uint8_t& _g, const uint8_t& _b) {
        uint64_t x = std::min((uint64_t)std::round(_x * scale_x), width-1);
        uint64_t t = _y*pix_per_path + link_pix_y;
        uint64_t s = t+pix_per_link;
        for (uint64_t y = t; y < s; ++y) {
            image[4 * width * y + 4 * x + 0] = _r;
            image[4 * width * y + 4 * x + 1] = _g;
            image[4 * width * y + 4 * x + 2] = _b;
            image[4 * width * y + 4 * x + 3] = 255;
        }
    };
    
    // map from path id to its starting y position
    //hash_map<uint64_t, uint64_t> path_layout_y;
    std::vector<uint64_t> path_layout_y; path_layout_y.resize(path_count);
    if (args::get(path_per_row)) {
        for (uint64_t i = 0; i < path_count; ++i) {
            path_layout_y[i] = i;
        }
    } else { // pack the paths
        // buffer to record layout bounds
        std::vector<bool> path_layout_buf; path_layout_buf.resize(path_count * width);
        std::vector<path_handle_t> path_order = algorithms::id_ordered_paths(graph);
        for (auto& path : path_order) {
            // get the block which this path covers
            uint64_t min_x = std::numeric_limits<uint64_t>::max();
            uint64_t max_x = std::numeric_limits<uint64_t>::min(); // 0
            graph.for_each_occurrence_in_path(path, [&](const occurrence_handle_t& occ) {
                    handle_t h = graph.get_occurrence(occ);
                    uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                    min_x = std::min(min_x, p);
                    max_x = std::max(max_x, p+graph.get_length(h));
                });
            //std::cerr << "min and max x " << min_x << " " << max_x << " vs " << len << std::endl;
            // now find where this would fit and mark the layout buffer
            // due to our sorted paths, we are able to drop to the lowest available layout position
            uint64_t path_y = 0;
            uint64_t actual_x = min_x*scale_x;
            while (path_layout_buf[width * path_y + actual_x] && path_y+1 < path_count) {
                ++path_y;
            }
            for (uint64_t i = min_x*scale_x; i < std::min((uint64_t)(max_x*scale_x + path_padding), width); ++i) {
                path_layout_buf[width * path_y + i] = 1;
            }
            //std::cerr << "path_y " << graph.get_path_name(path) << " " << path_count - path_y - 1 << std::endl;
            path_layout_y[as_integer(path)] = path_count - path_y - 1;
        }
    }

    graph.for_each_path_handle([&](const path_handle_t& path) {
            // use a sha256 to get a few bytes that we'll use for a color
            std::string path_name = graph.get_path_name(path);
            bool is_aln = false;
            if (aln_mode) {
                std::string::size_type n = path_name.find(aln_prefix);
                if (n == 0) {
                    is_aln = true;
                }
            }
            uint32_t path_name_hash = djb2_hash32(path_name.c_str());
            uint8_t path_r = 0;
            uint8_t path_b = 0;
            uint8_t path_g = 0;
            memcpy(&path_r, &path_name_hash, sizeof(uint8_t));
            memcpy(&path_b, ((uint8_t*)&path_name_hash)+1, sizeof(uint8_t));
            memcpy(&path_g, ((uint8_t*)&path_name_hash)+2, sizeof(uint8_t));
            float path_r_f = (float)path_r/(float)(std::numeric_limits<uint8_t>::max());
            float path_g_f = (float)path_g/(float)(std::numeric_limits<uint8_t>::max());
            float path_b_f = (float)path_b/(float)(std::numeric_limits<uint8_t>::max());
            float sum = path_r_f + path_g_f + path_b_f;
            path_r_f /= sum;
            path_g_f /= sum;
            path_b_f /= sum;
            if (is_aln) {
                float x = path_r_f;
                path_r_f = (x + 0.5*9)/10;
                path_g_f = (x + 0.5*9)/10;
                path_b_f = (x + 0.5*9)/10;
                // check the path orientations
                uint64_t steps = 0;
                uint64_t rev = 0;
                graph.for_each_occurrence_in_path(path, [&](const occurrence_handle_t& occ) {
                        handle_t h = graph.get_occurrence(occ);
                        ++steps;
                        rev += graph.get_is_reverse(h);
                    });
                bool is_rev = (float)rev/(float)steps > 0.5;
                if (is_rev) {
                    path_r_f = path_r_f * 0.9;
                    path_g_f = path_g_f * 0.9;
                    path_b_f = path_b_f * 1.2;
                } else {
                    path_b_f = path_b_f * 0.9;
                    path_g_f = path_g_f * 0.9;
                    path_r_f = path_r_f * 1.2;
                }
            }
            // brighten the color
            float f = std::min(1.5, 1.0/std::max(std::max(path_r_f, path_g_f), path_b_f));
            path_r = (uint8_t)std::round(255*std::min(path_r_f*f, (float)1.0));
            path_g = (uint8_t)std::round(255*std::min(path_g_f*f, (float)1.0));
            path_b = (uint8_t)std::round(255*std::min(path_b_f*f, (float)1.0));
            std::cerr << "path " << as_integer(path) << " " << graph.get_path_name(path) << " " << path_r_f << " " << path_g_f << " " << path_b_f
                      << " " << (int)path_r << " " << (int)path_g << " " << (int)path_b << std::endl;
            /// Loop over all the occurrences along a path, from first through last and draw them
            uint64_t path_rank = as_integer(path);
            uint64_t step = 0;
            graph.for_each_occurrence_in_path(path, [&](const occurrence_handle_t& occ) {
                    handle_t h = graph.get_occurrence(occ);
                    uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                    uint64_t hl = graph.get_length(h);
                    // make contects for the bases in the node
                    uint64_t path_y = path_layout_y[path_rank];
                    for (uint64_t i = 0; i < hl; i+=1/scale_x) {
                        add_path_step(p+i, path_y, path_r, path_g, path_b);
                    }
                 });
            // add in a visual motif that shows the links between path pieces
            // this is most meaningful in a linear layout
            if (args::get(link_path_pieces)) {
                uint64_t min_x = std::numeric_limits<uint64_t>::max();
                uint64_t max_x = std::numeric_limits<uint64_t>::min(); // 0
                graph.for_each_occurrence_in_path(path, [&](const occurrence_handle_t& occ) {
                        handle_t h = graph.get_occurrence(occ);
                        uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                        min_x = std::min(min_x, p);
                        max_x = std::max(max_x, p+graph.get_length(h));
                    });
                // now touch up the range
                uint64_t path_y = path_layout_y[path_rank];
                for (uint64_t i = min_x; i < max_x; i+=1/scale_x) {
                    add_path_link(i, path_y, path_r, path_g, path_b);
                }
            }
         });

    // trim vertical space to fit
    uint64_t min_y = std::numeric_limits<uint64_t>::max();
    uint64_t max_y = std::numeric_limits<uint64_t>::min(); // 0
    for (uint64_t y = 0; y < height + path_space; ++y) {
        for (uint64_t x = 0; x < width; ++x) {
            uint8_t r = image[4 * width * y + 4 * x + 0];
            uint8_t g = image[4 * width * y + 4 * x + 1];
            uint8_t b = image[4 * width * y + 4 * x + 2];
            uint8_t a = image[4 * width * y + 4 * x + 3];
            if (r != 255 || g != 255 || b != 255) {
                min_y = std::min(min_y, y);
                max_y = std::max(max_y, y);
            }
        }
    }
    // provide some default padding at the bottom, to clarify the edges
    max_y = std::min(path_space + height, max_y + bottom_padding);
    //std::cerr << "min and max y " << min_y << " " << max_y << std::endl;
    std::vector<uint8_t> crop;
    uint64_t crop_height = max_y - min_y;
    crop.resize(width * crop_height * 4, 255);
    for (uint64_t y = 0; y < max_y - min_y; ++y) {
        for (uint64_t x = 0; x < width; ++x) {
            crop[4 * width * y + 4 * x + 0] = image[4 * width * (y + min_y) + 4 * x + 0];
            crop[4 * width * y + 4 * x + 1] = image[4 * width * (y + min_y) + 4 * x + 1];
            crop[4 * width * y + 4 * x + 2] = image[4 * width * (y + min_y) + 4 * x + 2];
            crop[4 * width * y + 4 * x + 3] = image[4 * width * (y + min_y) + 4 * x + 3];
        }
    }
    
    png::encodeOneStep(filename, crop, width, crop_height);

    return 0;
}

static Subcommand odgi_viz("viz", "visualize the graph",
                           PIPELINE, 3, main_viz);


}
