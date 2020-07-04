#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "algorithms/hash.hpp"
#include "algorithms/id_ordered_paths.hpp"
#include "lodepng.h"
#include <limits>
#include "picosha2.h"
#include <iostream>

//#define debug_odgi_viz

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
        void encodeOneStep(const char *filename, std::vector<unsigned char> &image, unsigned width, unsigned height) {
            //Encode the image
            unsigned error = lodepng::encode(filename, image, width, height);

            //if there's an error, display it
            if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
        }

        //Example 2
        //Encode from raw pixels to an in-memory PNG file first, then write it to disk
        //The image argument has width * height RGBA pixels or width * height * 4 bytes
        void encodeTwoSteps(const char *filename, std::vector<unsigned char> &image, unsigned width, unsigned height) {
            std::vector<unsigned char> png;

            unsigned error = lodepng::encode(png, image, width, height);
            if (!error) lodepng::save_file(png, filename);

            //if there's an error, display it
            if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
        }

        //Example 3
        //Save a PNG file to disk using a State, normally needed for more advanced usage.
        //The image argument has width * height RGBA pixels or width * height * 4 bytes
        void encodeWithState(const char *filename, std::vector<unsigned char> &image, unsigned width, unsigned height) {
            std::vector<unsigned char> png;
            lodepng::State state; //optionally customize this one

            unsigned error = lodepng::encode(png, image, width, height, state);
            if (!error) lodepng::save_file(png, filename);

            //if there's an error, display it
            if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
        }

    }


    using namespace odgi::subcommand;

    int main_viz(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi viz";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser("variation graph visualizations");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the index from this file", {'i', "idx"});
        args::ValueFlag<std::string> png_out_file(parser, "FILE", "write the output (png) to this file", {'o', "out"});
        args::ValueFlag<uint64_t> image_width(parser, "N", "width in pixels of output image", {'x', "width"});
        args::ValueFlag<uint64_t> image_height(parser, "N", "height in pixels of output image", {'y', "height"});
        args::ValueFlag<uint64_t> path_height(parser, "N", "path display height", {'P', "path-height"});
        args::ValueFlag<uint64_t> path_x_pad(parser, "N", "path x padding", {'X', "path-x-padding"});
        args::Flag pack_paths(parser, "bool", "pack the graphs rather than displaying a single path per row",{'R', "pack-paths"});
        args::ValueFlag<float> link_path_pieces(parser, "FLOAT","show thin links of this relative width to connect path pieces",{'L', "link-path-pieces"});
        args::ValueFlag<std::string> alignment_prefix(parser, "STRING","apply alignment-related visual motifs to paths with this name prefix",{'A', "alignment-prefix"});
        args::Flag show_strands(parser, "bool","use reds and blues to show forward and reverse alignments (depends on -A)",{'S', "show-strand"});
        args::Flag binned_mode(parser, "binned-mode", "bin the variation graph before its visualization", {'b', "binned-mode"});
        args::ValueFlag<uint64_t> bin_width(parser, "bp", "width of each bin in basepairs along the graph vector",{'w', "bin-width"});
        args::Flag drop_gap_links(parser, "drop-gap-links", "don't include gap links in the output", {'g', "no-gap-links"});
        args::Flag color_by_nt_pos(parser, "color-by-nt-pos", "change the color intensity based on nucleotide positionb", {'c', "color-by-nt-pos"});
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
        if (argc == 1) {
            std::cout << parser;
            return 1;
        }

        if (args::get(threads)) {
            omp_set_num_threads(args::get(threads));
        } else {
            omp_set_num_threads(1);
        }

        if (!dg_in_file) {
            std::cerr
                    << "[odgi viz] error: Please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
                    << std::endl;
            return 1;
        }

        if (!png_out_file) {
            std::cerr
                    << "[odgi viz] error: Please specify an output file to where to store the PNG via -o=[FILE], --out=[FILE]."
                    << std::endl;
            return 1;
        }

        if (!args::get(binned_mode) && ((args::get(bin_width) > 0) || (args::get(drop_gap_links)))){
            std::cerr
                    << "[odgi viz] error: Please specify the -b/--binned-mode option to use the -w/--bin_width and -g/--no-gap-links options."
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

        //NOTE: this sample will overwrite the file or test.png without warning!
        //const char* filename = argc > 1 ? argv[1] : "test.png";
        if (!args::get(png_out_file).size()) {
            std::cerr << "[odgi viz] error: output image required" << std::endl;
            return 1;
        }
        const char *filename = args::get(png_out_file).c_str();

        std::vector<uint64_t> position_map(graph.get_node_count() + 1);
        uint64_t len = 0;
        nid_t last_node_id = graph.min_node_id();
        graph.for_each_handle([&](const handle_t &h) {
            nid_t node_id = graph.get_id(h);
            if (node_id - last_node_id > 1) {
                std::cerr << "[odgi viz] error: The graph is not optimized. Please run 'odgi sort' using -O, --optimize" << std::endl;
                exit(1);
            }
            last_node_id = node_id;

            position_map[number_bool_packing::unpack_number(h)] = len;
            uint64_t hl = graph.get_length(h);
            len += hl;

#ifdef debug_odgi_viz
            std::cerr << "SEGMENT ID: " << graph.get_id(h) << " - " << as_integer(h) << " - index_in_position_map (" << number_bool_packing::unpack_number(h) << ") = " << len << std::endl;
#endif
        });
        position_map[position_map.size() - 1] = len;

        uint64_t path_count = graph.get_path_count();
        uint64_t pix_per_path = (args::get(path_height) ? args::get(path_height) : 10);
        uint64_t pix_per_link = std::max((uint64_t) 1, (uint64_t) std::round(args::get(link_path_pieces) * pix_per_path));
        uint64_t link_pix_y = pix_per_path / 2 - pix_per_link / 2;
        uint64_t path_space = path_count * pix_per_path;
        uint64_t path_padding = args::get(path_x_pad);
        uint64_t bottom_padding = 5;
        // the math here only works if the image size is greater than or equal to the graph length
        uint64_t width = std::min(len, (args::get(image_width) ? args::get(image_width) : 1000));
        uint64_t height = std::min(len, (args::get(image_height) ? args::get(image_height) : 1000) + bottom_padding);
        float scale_x = (float) width / (float) len;
        float scale_y = (float) height / (float) len;

        float _bin_width = args::get(bin_width);
        if (args::get(binned_mode)){
            if (_bin_width == 0){
                _bin_width = 1 / scale_x; // It can be a float value.
            }else{
                float num_bins = (float) len / _bin_width;// + (len % bin_width ? 1 : 0);
                width = num_bins;
            }

            // Each pixel corresponds to a bin
            scale_x = 1; //scale_x*bin_width;

            std::cerr << "Binned mode" << std::endl;
            std::cerr << "bin width: " << _bin_width << std::endl;
            std::cerr << "image width: " << width << std::endl;
        }else{
            _bin_width = 1;
        }

        if (width > 50000){
            std::cerr
                    << "[odgi viz] warning: you are going to create a big image (width > 50000 pixels)."
                    << std::endl;
        }

        std::vector<uint8_t> image;
        image.resize(width * (height + path_space) * 4, 255);

        bool aln_mode = !args::get(alignment_prefix).empty();
        std::string aln_prefix;
        if (aln_mode) {
            aln_prefix = args::get(alignment_prefix);
        }

        auto add_point = [&](const uint64_t &_x, const uint64_t &_y,
                             const uint8_t &_r, const uint8_t &_g, const uint8_t &_b) {
            uint64_t x = std::min((uint64_t) std::round(_x * scale_x), width - 1);
            uint64_t y = std::min((uint64_t) std::round(_y * scale_y), height - 1) + path_space;
            image[4 * width * y + 4 * x + 0] = _r;
            image[4 * width * y + 4 * x + 1] = _g;
            image[4 * width * y + 4 * x + 2] = _b;
            image[4 * width * y + 4 * x + 3] = 255;
        };

        auto add_edge_from_positions = [&](uint64_t a, const uint64_t b) {
#ifdef debug_odgi_viz
            std::cerr << "Edge displayed" << std::endl;
            std::cerr << a << " --> " << b << std::endl;
#endif
            // In binned mode, the Links have to be tall to be visible; in standard mode, _bin_width is 1, so nothing changes here
            uint64_t dist = (b - a)*_bin_width;

            uint64_t i = 0;

            for (; i < dist; i += 1 / scale_y) {
                add_point(a, i, 0, 0, 0);
            }
            while (a <= b) {
                add_point(a, i, 0, 0, 0);
                a += 1 / scale_x;
            }
            for (uint64_t j = 0; j < dist; j += 1 / scale_y) {
                add_point(b, j, 0, 0, 0);
            }
        };

        auto add_edge_from_handles = [&](const handle_t& h, const handle_t& o) {
            uint64_t _a = position_map[number_bool_packing::unpack_number(h) + !number_bool_packing::unpack_bit(h)];
            uint64_t _b = position_map[number_bool_packing::unpack_number(o) + number_bool_packing::unpack_bit(o)];
            uint64_t a = std::min(_a, _b);
            uint64_t b = std::max(_a, _b);

#ifdef debug_odgi_viz
            std::cerr << graph.get_id(h) << " (" << number_bool_packing::unpack_bit(h) << ") --> " << b << " (" << number_bool_packing::unpack_bit(h) << ") " << std::endl;
#endif

            add_edge_from_positions(a, b);
        };

        if (args::get(binned_mode)){
            graph.for_each_handle([&](const handle_t &h) {
                uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                uint64_t hl = graph.get_length(h);

                int64_t last_bin = 0; // flag meaning "null bin"
                // make contents for the bases in the node
                for (uint64_t k = 0; k < hl; ++k) {
                    int64_t curr_bin = (p + k) / _bin_width + 1;
                    if (curr_bin != last_bin) {
#ifdef debug_odgi_viz
                        std::cerr << "position in map (" << p  << ") - curr_bin: " << curr_bin << std::endl;
#endif
                        add_point(curr_bin - 1, 0, 0, 0, 0);
                    }

                    last_bin = curr_bin;
                }
            });

            // The links are created later after the binning step.
        }else{
            graph.for_each_handle([&](const handle_t &h) {
                uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                uint64_t hl = graph.get_length(h);
                // make contents for the bases in the node
                //for (uint64_t i = 0; i < hl; ++i) {
                for (uint64_t i = 0; i < hl; i += 1 / scale_x) {
                    add_point(p + i, 0, 0, 0, 0);
                }
            });

            graph.for_each_handle([&](const handle_t& h) {
                // add contacts for the edges
                graph.follow_edges(h, false, [&](const handle_t& o) {
                    add_edge_from_handles(h, o);
                });
                graph.follow_edges(h, true, [&](const handle_t& o) {
                    add_edge_from_handles(o, h);
                });
            });
        }


        auto add_path_step = [&](const uint64_t &_x, const uint64_t &_y,
                                 const uint8_t &_r, const uint8_t &_g, const uint8_t &_b) {
            uint64_t x = std::min((uint64_t) std::round(_x * scale_x), width - 1);
            uint64_t t = _y * pix_per_path;
            uint64_t s = t + pix_per_path;
            for (uint64_t y = t; y < s; ++y) {
                image[4 * width * y + 4 * x + 0] = _r;
                image[4 * width * y + 4 * x + 1] = _g;
                image[4 * width * y + 4 * x + 2] = _b;
                image[4 * width * y + 4 * x + 3] = 255;
            }
        };

        auto add_path_link = [&](const uint64_t &_x, const uint64_t &_y,
                                 const uint8_t &_r, const uint8_t &_g, const uint8_t &_b) {
            uint64_t x = std::min((uint64_t) std::round(_x * scale_x), width - 1);
            uint64_t t = _y * pix_per_path + link_pix_y;
            uint64_t s = t + pix_per_link;
            for (uint64_t y = t; y < s; ++y) {
                image[4 * width * y + 4 * x + 0] = _r;
                image[4 * width * y + 4 * x + 1] = _g;
                image[4 * width * y + 4 * x + 2] = _b;
                image[4 * width * y + 4 * x + 3] = 255;
            }
        };

        // map from path id to its starting y position
        //hash_map<uint64_t, uint64_t> path_layout_y;
        std::vector<uint64_t> path_layout_y;
        path_layout_y.resize(path_count);
        if (!args::get(pack_paths)) {
            for (uint64_t i = 0; i < path_count; ++i) {
                path_layout_y[i] = i;
            }
        } else { // pack the paths
            // buffer to record layout bounds
            std::vector<bool> path_layout_buf;
            path_layout_buf.resize(path_count * width);
            std::vector<path_handle_t> path_order = algorithms::id_ordered_paths(graph);
            for (auto &path : path_order) {
                // get the block which this path covers
                uint64_t min_x = std::numeric_limits<uint64_t>::max();
                uint64_t max_x = std::numeric_limits<uint64_t>::min(); // 0
                graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                    handle_t h = graph.get_handle_of_step(occ);
                    uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                    min_x = std::min(min_x, p);
                    max_x = std::max(max_x, p + graph.get_length(h));
                });
                //std::cerr << "min and max x " << min_x << " " << max_x << " vs " << len << std::endl;
                // now find where this would fit and mark the layout buffer
                // due to our sorted paths, we are able to drop to the lowest available layout position
                uint64_t path_y = 0;
                uint64_t actual_x = min_x * scale_x;
                while (path_layout_buf[width * path_y + actual_x] && path_y + 1 < path_count) {
                    ++path_y;
                }
                for (uint64_t i = min_x * scale_x;
                     i < std::min((uint64_t) (max_x * scale_x + path_padding), width); ++i) {
                    path_layout_buf[width * path_y + i] = 1;
                }
                //std::cerr << "path_y " << graph.get_path_name(path) << " " << path_count - path_y - 1 << std::endl;
                path_layout_y[as_integer(path) - 1] = path_count - path_y - 1;
            }
        }

        std::unordered_set<pair<uint64_t, uint64_t>> edges_drawn;
        uint64_t gap_links_removed = 0;
        uint64_t total_links = 0;
        graph.for_each_path_handle([&](const path_handle_t &path) {
            // use a sha256 to get a few bytes that we'll use for a color
            std::string path_name = graph.get_path_name(path);

#ifdef debug_odgi_viz
            std::cerr << "path_name: " << path_name << std::endl;
#endif

            bool is_aln = false;
            if (aln_mode) {
                std::string::size_type n = path_name.find(aln_prefix);
                if (n == 0) {
                    is_aln = true;
                }
            }
            // use a sha256 to get a few bytes that we'll use for a color
            picosha2::byte_t hashed[picosha2::k_digest_size];
            picosha2::hash256(path_name.begin(), path_name.end(), hashed, hashed + picosha2::k_digest_size);
            uint8_t path_r = hashed[24];
            uint8_t path_g = hashed[8];
            uint8_t path_b = hashed[16];
            float path_r_f = (float) path_r / (float) (std::numeric_limits<uint8_t>::max());
            float path_g_f = (float) path_g / (float) (std::numeric_limits<uint8_t>::max());
            float path_b_f = (float) path_b / (float) (std::numeric_limits<uint8_t>::max());
            float sum = path_r_f + path_g_f + path_b_f;
            path_r_f /= sum;
            path_g_f /= sum;
            path_b_f /= sum;

            bool _color_by_nt_pos = args::get(color_by_nt_pos);

            // Calculate the number or steps, the reverse steps and the length of the path if any of this information
            // is needed depending on the input arguments.
            uint64_t steps = 0;
            uint64_t rev = 0;
            uint64_t path_len = 0;
            if ((is_aln && args::get(show_strands)) || _color_by_nt_pos){
                graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                    handle_t h = graph.get_handle_of_step(occ);
                    ++steps;
                    rev += graph.get_is_reverse(h);

                    path_len += graph.get_length(h);
                });
            }

            if (is_aln) {
                float x = path_r_f;
                path_r_f = (x + 0.5 * 9) / 10;
                path_g_f = (x + 0.5 * 9) / 10;
                path_b_f = (x + 0.5 * 9) / 10;
                // check the path orientations
                if (args::get(show_strands)) {
                    bool is_rev = (float) rev / (float) steps > 0.5;
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
            }

            // brighten the color
            float f = std::min(1.5, 1.0 / std::max(std::max(path_r_f, path_g_f), path_b_f));
            path_r = (uint8_t) std::round(255 * std::min(path_r_f * f, (float) 1.0));
            path_g = (uint8_t) std::round(255 * std::min(path_g_f * f, (float) 1.0));
            path_b = (uint8_t) std::round(255 * std::min(path_b_f * f, (float) 1.0));
            //std::cerr << "path " << as_integer(path) << " " << graph.get_path_name(path) << " " << path_r_f << " " << path_g_f << " " << path_b_f
            //          << " " << (int)path_r << " " << (int)path_g << " " << (int)path_b << std::endl;
            uint64_t path_rank = as_integer(path) - 1;

            uint64_t curr_len = 0;
            float x = 1;
            if (args::get(binned_mode)) {
                std::vector<std::pair<uint64_t, uint64_t>> links;
                std::vector<uint64_t> bin_ids;
                int64_t last_bin = 0; // flag meaning "null bin"
                graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                    handle_t h = graph.get_handle_of_step(occ);
                    uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                    uint64_t hl = graph.get_length(h);

                    // make contents for the bases in the node

                    uint64_t path_y = path_layout_y[path_rank];
                    for (uint64_t k = 0; k < hl; ++k) {
                        int64_t curr_bin = (p + k) / _bin_width + 1;
                        if (curr_bin != last_bin) {
#ifdef debug_odgi_viz
                            std::cerr << "curr_bin: " << curr_bin << std::endl;
#endif

                            if (_color_by_nt_pos){
                                x = 1 - ( (float)(curr_len + k) / (float)(len))*0.8;
                            }
                            add_path_step(curr_bin - 1, path_y, (float)path_r * x, (float)path_g * x, (float)path_b * x);

                            if (std::abs(curr_bin - last_bin) > 1 || last_bin == 0) {
                                // bin cross!
                                links.push_back(std::make_pair(last_bin, curr_bin));
                            }
                        }
                        bin_ids.push_back(curr_bin);

                        last_bin = curr_bin;
                    }

                    curr_len += hl;
                });
                links.push_back(std::make_pair(last_bin, 0));

                if (args::get(drop_gap_links)) {
                    std::sort(bin_ids.begin(), bin_ids.end());
                    total_links += links.size();

                    uint64_t fill_pos = 0;

                    for (uint64_t i = 0; i < links.size(); ++i) {
                        auto link = links[i];

                        if (link.first == 0 || link.second == 0)
                            continue;

                        if (link.first > link.second) {
                            links[fill_pos++] = link;
                            continue;
                        }

                        auto left_it = std::lower_bound(bin_ids.begin(), bin_ids.end(), link.first + 1);
                        auto right_it = std::lower_bound(bin_ids.begin(), bin_ids.end(), link.second);
                        if (right_it > left_it) {
                            links[fill_pos++] = link;
                        }
                    }

                    gap_links_removed += links.size() - fill_pos;
                    links.resize(fill_pos);
                }

                for (auto const link : links) {
#ifdef debug_odgi_viz
                    std::cerr << link.first << " --> " << link.second << std::endl;
#endif
                    if( (link.first != 0) && (link.second != 0) && std::abs((int64_t)link.first - (int64_t)link.second) > 1){
                        uint64_t a = std::min(link.first, link.second) - 1;
                        uint64_t b = std::max(link.first, link.second) - 1;

                        auto pair = make_pair(a, b);

                        // Check if the edge is already displayed
                        if (edges_drawn.find(pair) == edges_drawn.end()) {
                            edges_drawn.insert(pair);

                            // add contents for the edge
                            add_edge_from_positions(a, b);
                        }
                    }
                }
            }else{
                /// Loop over all the steps along a path, from first through last and draw them
                graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
                    handle_t h = graph.get_handle_of_step(occ);
                    uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                    uint64_t hl = graph.get_length(h);
                    // make contects for the bases in the node
                    uint64_t path_y = path_layout_y[path_rank];
                    for (uint64_t i = 0; i < hl; i+=1/scale_x) {
                        if (_color_by_nt_pos){
                            x = 1 - ((float)(curr_len + i*scale_x) / (float)(len))*0.8;
                        }
                        add_path_step(p+i, path_y, (float)path_r * x, (float)path_g * x, (float)path_b * x);
                    }

                    curr_len += hl;
                });
            }

            // add in a visual motif that shows the links between path pieces
            // this is most meaningful in a linear layout
            if (args::get(link_path_pieces)) {
                uint64_t min_x = std::numeric_limits<uint64_t>::max();
                uint64_t max_x = std::numeric_limits<uint64_t>::min(); // 0

                // In binned mode, the min/max_x values changes based on the bin width; in standard mode, _bin_width is 1, so nothing changes here
                graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                    handle_t h = graph.get_handle_of_step(occ);
                    uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                    min_x = std::min(min_x, (uint64_t)(p / _bin_width));
                    max_x = std::max(max_x, (uint64_t)((p + graph.get_length(h)) / _bin_width));
                });

                // now touch up the range
                uint64_t path_y = path_layout_y[path_rank];
                for (uint64_t i = min_x; i < max_x; i += 1 / scale_x) {
                    add_path_link(i, path_y, path_r, path_g, path_b);
                }
            }
        });

        if (args::get(drop_gap_links)) {
            std::cerr << "Gap links removed: " << gap_links_removed << " of " << total_links << " total links" << std::endl;
        }

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

    static Subcommand odgi_viz("viz", "visualize the graph",PIPELINE, 3, main_viz);

}