#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/bin_path_info.hpp"
#include "algorithms/hash.hpp"
#include "algorithms/id_ordered_paths.hpp"
#include "lodepng.h"
#include <limits>
#include <regex>
#include "picosha2.h"
#include "algorithms/draw.hpp"

//#define debug_odgi_viz

#define RGB_BIN_LINKS 110

#include "fonts/font5x8.h"


#define PATH_NAMES_MAX_NUM_OF_CHARACTERS 128
#define PATH_NAMES_MAX_CHARACTER_SIZE 64


namespace odgi {

    using namespace odgi::subcommand;

    // helper to get the prefix of a string
    std::string prefix(const std::string& s, const char c) {
        //std::cerr << "prefix of " << s << " by " << c << " is " << s.substr(0, s.find(c)) << std::endl;
        return s.substr(0, s.find(c));
    }

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
        args::ValueFlag<std::string> alignment_prefix(parser, "STRING","apply alignment-related visual motifs to paths with this name prefix (it affects the -S and -d options)",{'A', "alignment-prefix"});
        args::Flag show_strands(parser, "bool","use reds and blues to show forward and reverse alignments",{'S', "show-strand"});
        args::Flag color_by_mean_inversion_rate(parser, "color-by-mean-inversion-rate", "change the color respect to the node strandness (black for forward, red for reverse); in binned mode, change the color respect to the mean inversion rate of the path for each bin, from black (no inversions) to red (bin mean inversion rate equals to 1)", {'z', "color-by-mean-inversion-rate"});

        args::ValueFlag<char> _color_by_prefix(parser, "C", "colors paths by their names looking at the prefix before the given character C",{'s', "color-by-prefix"});

        /// Range selection
        args::ValueFlag<std::string> _nucleotide_pangenomic_range(parser, "STRING","nucleotide pangenomic range to visualize [pos1-pos2]",{'r', "pangenomic-range"});

        /// Paths selection
        args::ValueFlag<std::string> path_names_file(parser, "FILE", "list of paths to display in the specified order; the file must contain one path name per line and a subset of all paths can be specified.", {'p', "paths-to-display"});

        /// Path names
        args::Flag hide_path_names(parser, "bool", "hide path names on the left",{'H', "hide-path-names"});
        args::Flag color_path_names_background(parser, "bool", "color path names background with the same color as paths",{'C', "color-path-names-background"});
        args::ValueFlag<uint64_t> _max_num_of_characters(parser, "N", "max number of characters to display for each path name (default: 20)",{'c', "max-num-of-characters"});

        /// Binned mode
        args::Flag binned_mode(parser, "binned-mode", "bin the variation graph before its visualization", {'b', "binned-mode"});
        args::ValueFlag<uint64_t> bin_width(parser, "bp", "width of each bin in basepairs along the graph vector",{'w', "bin-width"});
        args::Flag drop_gap_links(parser, "drop-gap-links", "don't include gap links in the output", {'g', "no-gap-links"});
        args::Flag color_by_mean_coverage(parser, "color-by-mean-coverage", "change the color respect to the mean coverage of the path for each bin, from black (no coverage) to blue (max bin mean coverage in the entire graph)", {'m', "color-by-mean-coverage"});

        /// Gradient mode
        args::Flag change_darkness(parser, "change-darkness", "change the color darkness based on nucleotide position in the path", {'d', "change-darkness"});
        args::Flag longest_path(parser, "longest-path", "use the longest path length to change the color darkness", {'l', "longest-path"});
        args::Flag white_to_black(parser, "white-to-black", "change the color darkness from white (for the first nucleotide position) to black (for the last nucleotide position)", {'u', "white-to-black"});

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

        if (
                !args::get(binned_mode) &&
                ((args::get(bin_width) > 0) || args::get(drop_gap_links) ||
                args::get(color_by_mean_coverage))
                ){
            std::cerr
                    << "[odgi viz] error: Please specify the -b/--binned-mode option to use the "
                       "-w/--bin_width, -g/--no-gap-links, and -m/--color-by-mean-coverage "
                       "options."
                    << std::endl;
            return 1;
        }

        if (!args::get(change_darkness) && (args::get(longest_path) || args::get(white_to_black))){
            std::cerr
                    << "[odgi viz] error: Please specify the -d/--change-darkness option to use the -l/--longest-path and -u/--white-to-black options."
                    << std::endl;
            return 1;
        }

        if ((args::get(_color_by_prefix) != 0) + args::get(show_strands) + args::get(white_to_black) + args::get(color_by_mean_coverage) + args::get(color_by_mean_inversion_rate) > 1) {
            std::cerr
                    << "[odgi viz] error: Please specify only one of the following options: "
                       "-s/--color-by-prefix, -S/--show-strand, -u/--white-to-black, "
                       "-m/--color-by-mean-coverage, and -z/--color-by-mean-inversion."
                    << std::endl;
            return 1;
        }

        if (args::get(change_darkness) && (args::get(color_by_mean_coverage) || args::get(color_by_mean_inversion_rate))) {
            std::cerr
                    << "[odgi viz] error: Please specify the -d/--change-darkness option without specifying "
                       "-m/--color-by-mean-coverage or -z/--color-by-mean-inversion."
                    << std::endl;
            return 1;
        }

        if (args::get(pack_paths) && !args::get(path_names_file).empty()){
            std::cerr
                    << "[odgi viz] error: Please specify -R/--pack-paths or -p/--paths-to-display, not both."
                    << std::endl;
            return 1;
        }

        if (args::get(hide_path_names) && (args::get(color_path_names_background) || args::get(_max_num_of_characters))){
            std::cerr
                    << "[odgi viz] error: Please specify the -C/--color-path-names-background and -c/--max-num-of-characters "
                       "options without specifying -H/--hide-path-names." << std::endl;
            return 1;
        }


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

        //NOTE: this sample will overwrite the file or test.png without warning!
        //const char* filename = argc > 1 ? argv[1] : "test.png";
        if (args::get(png_out_file).empty()) {
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

        double pangenomic_start_pos = 0.0;
        double pangenomic_end_pos = (double) (len - 1);

        {
            std::string nucleotide_pangenomic_range = args::get(_nucleotide_pangenomic_range);
            if (!nucleotide_pangenomic_range.empty()) {
                std::regex regex("-");
                std::vector<std::string> splitted(
                        std::sregex_token_iterator(args::get(_nucleotide_pangenomic_range).begin(), args::get(_nucleotide_pangenomic_range).end(), regex, -1),
                        std::sregex_token_iterator()
                );


                if (splitted.size() != 2) {
                    //todo to improve the message
                    std::cerr
                            << "[odgi viz] error: Please specify a nucleotide pangenomic range with a start position "
                               "and an end position."
                            << std::endl;
                    return 1;
                }

                //todo checks that they are numbers

                pangenomic_start_pos = std::min((double)len - 1, (double)stoi(splitted[0]));
                pangenomic_end_pos = std::min((double)len - 1, (double)stoi(splitted[1]));

                if (pangenomic_start_pos >= pangenomic_end_pos) {
                    std::cerr
                            << "[odgi viz] error: Please specify a start position less than the end position."
                            << std::endl;
                    return 1;
                }
            }
        }

        uint64_t max_num_of_characters = args::get(_max_num_of_characters) > 1 ? min(args::get(_max_num_of_characters), (uint64_t) PATH_NAMES_MAX_NUM_OF_CHARACTERS) : 32;
        uint64_t path_count = graph.get_path_count();
        uint64_t pix_per_path = args::get(path_height) ? args::get(path_height) : 10;
        uint64_t pix_per_link = std::max((uint64_t) 1, (uint64_t) std::round(args::get(link_path_pieces) * pix_per_path));
        uint64_t link_pix_y = pix_per_path / 2 - pix_per_link / 2;
        uint64_t path_padding = args::get(path_x_pad);
        uint64_t bottom_padding = 5;

        uint64_t len_to_visualize = pangenomic_end_pos - pangenomic_start_pos + 1;

        // the math here only works if the image size is greater than or equal to the graph length
        uint64_t width = std::min(len_to_visualize, (args::get(image_width) ? args::get(image_width) : 1000));
        uint64_t height = std::min(len_to_visualize, (args::get(image_height) ? args::get(image_height) : 1000) + bottom_padding);

        float scale_x = (float) width / (float) len_to_visualize;
        float scale_y = (float) height / (float) len_to_visualize;

        float _bin_width = args::get(bin_width);
        bool _binned_mode = args::get(binned_mode);
        if (_binned_mode){
            if (_bin_width == 0){
                _bin_width = 1 / scale_x; // It can be a float value.
            }else{
                width = len_to_visualize / _bin_width;// + (len_to_visualize % bin_width ? 1 : 0);
            }

            pangenomic_start_pos /= _bin_width;
            pangenomic_end_pos /= _bin_width;

            // Each pixel corresponds to a bin
            scale_x = 1; //scale_x*bin_width;

            std::cerr << "Binned mode" << std::endl;
            std::cerr << "bin width: " << _bin_width << std::endl;
            std::cerr << "image width: " << width << std::endl;
        }else{
            _bin_width = 1;
        }

        /*std::cerr << "real len: " << len << std::endl;
        std::cerr << "pangenomic_start_pos: " << pangenomic_start_pos << "\npangenomic_end_pos: " << pangenomic_end_pos << std::endl;
        std::cerr << "len_to_visualize: " << len_to_visualize << std::endl;*/

        // After the bin_width and scale_xy calculations
        uint16_t width_path_names = 0;
        uint8_t max_num_of_chars = 0;
        uint8_t char_size = 0;

        // map from path id to its starting y position
        //hash_map<uint64_t, uint64_t> path_layout_y;
        std::vector<int64_t> path_layout_y;
        path_layout_y.resize(path_count, -1);
        if (!args::get(pack_paths)) {
            std::string _path_names = args::get(path_names_file);
            if (!_path_names.empty()){
                std::ifstream path_names_in(_path_names);

                uint64_t rank_for_visualization = 0;
                uint64_t num_of_paths_in_file = 0;

                std::string line;
                while (std::getline(path_names_in, line)) {
                    if (!line.empty()){
                        if (graph.has_path(line)){
                            uint64_t path_rank = as_integer(graph.get_path_handle(line)) - 1;
                            if (path_layout_y[path_rank] < 0){
                                path_layout_y[path_rank] = rank_for_visualization;
                            }else{
                                std::cerr << "[odgi viz] error: In the path list there are duplicated path names." << std::endl;
                                exit(1);
                            }

                            rank_for_visualization++;
                        }

                        num_of_paths_in_file++;
                    }
                }

                std::cerr << "Found " << rank_for_visualization << "/" << num_of_paths_in_file << " paths to display." << std::endl;

                if (rank_for_visualization == 0){
                    std::cerr << "[odgi viz] error: No path to display." << std::endl;
                    exit(1);
                }

                path_count = rank_for_visualization;
            }else{
                for (uint64_t i = 0; i < path_count; ++i) {
                    path_layout_y[i] = i;
                }
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
                //std::cerr << "min and max x " << min_x << " " << max_x << " vs " << len_to_visualize << std::endl;
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

        uint64_t path_space = path_count * pix_per_path;

        std::vector<uint8_t> image;
        image.resize(width * (height + path_space) * 4, 255);

        std::vector<uint8_t> image_path_names;
        if (!args::get(hide_path_names) && !args::get(pack_paths) && pix_per_path >= 8) {
            uint64_t _max_num_of_chars = std::numeric_limits<uint64_t>::min();
            graph.for_each_path_handle([&](const path_handle_t &path) {
                uint64_t path_rank = as_integer(path) - 1;
                if (path_layout_y[path_rank] >= 0){
                    _max_num_of_chars = max((uint64_t) _max_num_of_chars, graph.get_path_name(path).length());
                }
            });
            max_num_of_chars = min(_max_num_of_chars, max_num_of_characters);

            char_size = min((uint16_t)((pix_per_path / 8) * 8), (uint16_t) PATH_NAMES_MAX_CHARACTER_SIZE);

            width_path_names = max_num_of_chars * char_size + char_size / 2;

            image_path_names.resize(width_path_names * (height + path_space) * 4, 255);
        }

        if (width_path_names + width > 50000){
            std::cerr
                    << "[odgi viz] warning: you are going to create a big image (width > 50000 pixels)."
                    << std::endl;
        }

        bool aln_mode = !args::get(alignment_prefix).empty();
        std::string aln_prefix;
        if (aln_mode) {
            aln_prefix = args::get(alignment_prefix);
        }

        char path_name_prefix_separator = '\0';
        if (_color_by_prefix) {
            path_name_prefix_separator = args::get(_color_by_prefix);
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

        auto add_edge_from_positions = [&](uint64_t a, const uint64_t b, uint8_t rgb) {
#ifdef debug_odgi_viz
            std::cerr << "Edge displayed" << std::endl;
            std::cerr << a << " --> " << b << std::endl;
#endif
            // In binned mode, the Links have to be tall to be visible; in standard mode, _bin_width is 1, so nothing changes here
            uint64_t dist = (b - a) * _bin_width;

            uint64_t i = 0;

            for (; i < dist; i += 1 / scale_y) {
                if (a >= pangenomic_start_pos && a <= pangenomic_end_pos) {
                    add_point(a - pangenomic_start_pos, i, rgb, rgb, rgb);
                }
            }

            while (a <= b) {
                if (a >= pangenomic_start_pos && a <= pangenomic_end_pos) {
                    add_point(a - pangenomic_start_pos, i, rgb, rgb, rgb);
                }
                a += 1 / scale_x;
            }
            if (b >= pangenomic_start_pos && b <= pangenomic_end_pos) {
                for (uint64_t j = 0; j < dist; j += 1 / scale_y) {
                    add_point(b - pangenomic_start_pos, j, rgb, rgb, rgb);
                }
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

            add_edge_from_positions(a, b, 0);
        };

        {
            std::function<bool(const handle_t)> is_a_handle_to_hide;
            if (path_count < graph.get_path_count()){
                is_a_handle_to_hide = [&](const handle_t &h) {
                    return graph.for_each_step_on_handle(h, [&](const step_handle_t &step) {
                        auto path_handle = graph.get_path_handle_of_step(step);
                        if (path_layout_y[as_integer(path_handle) - 1] >= 0) {
                            return false;
                        }

                        return true;
                    });
                };
            }else{
                is_a_handle_to_hide = [&](const handle_t &h) {
                    return false;
                };
            }

            if (_binned_mode){
                graph.for_each_handle([&](const handle_t &h) {
                    if (!is_a_handle_to_hide(h)){
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
                                if (curr_bin - 1 >= pangenomic_start_pos && curr_bin - 1 <= pangenomic_end_pos) {
                                    add_point(curr_bin - 1 - pangenomic_start_pos, 0, RGB_BIN_LINKS, RGB_BIN_LINKS, RGB_BIN_LINKS);
                                }
                            }

                            last_bin = curr_bin;
                        }
                    }
                });

                // The links are created later after the binning step.
            }else{
                graph.for_each_handle([&](const handle_t &h) {
                    if (!is_a_handle_to_hide(h)){
                        uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                        uint64_t hl = graph.get_length(h);
                        // make contents for the bases in the node
                        for (uint64_t i = 0; i < hl; i += 1 / scale_x) {
                            if ((p + i) >= pangenomic_start_pos && (p + i) <= pangenomic_end_pos) {
                                add_point(p + i - pangenomic_start_pos, 0, 0, 0, 0);
                            }
                        }

                        // add contacts for the edges
                        graph.follow_edges(h, false, [&](const handle_t& o) {
                            if (!is_a_handle_to_hide(o)){
                                add_edge_from_handles(h, o);
                            }
                        });
                        graph.follow_edges(h, true, [&](const handle_t& o) {
                            if (!is_a_handle_to_hide(o)){
                                add_edge_from_handles(o, h);
                            }
                        });
                    }
                });
            }
        }

        auto add_path_step = [&](std::vector<uint8_t> &img, uint64_t width_img,
                const uint64_t &_x, const uint64_t &_y,
                const uint8_t &_r, const uint8_t &_g, const uint8_t &_b) {
            uint64_t x = std::min((uint64_t) std::round(_x * scale_x), width - 1);
            uint64_t t = _y * pix_per_path;
            uint64_t s = t + pix_per_path;

            for (uint64_t y = t; y < s; ++y) {
                img[4 * width_img * y + 4 * x + 0] = _r;
                img[4 * width_img * y + 4 * x + 1] = _g;
                img[4 * width_img * y + 4 * x + 2] = _b;
                img[4 * width_img * y + 4 * x + 3] = 255;
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

        bool _show_strands = args::get(show_strands);

        bool _change_darkness = args::get(change_darkness);
        bool _longest_path = args::get(longest_path);
        bool _white_to_black = args::get(white_to_black);

        bool _color_by_mean_coverage = args::get(color_by_mean_coverage);
        bool _color_by_mean_inversion_rate = args::get(color_by_mean_inversion_rate);

        uint64_t longest_path_len = 0;
        double max_mean_cov = 0.0;
        if ((_change_darkness && _longest_path) || (_binned_mode && _color_by_mean_coverage)){
            graph.for_each_path_handle([&](const path_handle_t &path) {
                uint64_t path_rank = as_integer(path) - 1;
                if (path_layout_y[path_rank] >= 0){
                    uint64_t curr_len = 0, p, hl;
                    handle_t h;
                    std::map<uint64_t, algorithms::path_info_t> bins;
                    graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                        h = graph.get_handle_of_step(occ);
                        hl = graph.get_length(h);

                        curr_len += hl;

                        p = position_map[number_bool_packing::unpack_number(h)];
                        for (uint64_t k = 0; k < hl; ++k) {
                            int64_t curr_bin = (p + k) / _bin_width + 1;

                            ++bins[curr_bin].mean_cov;
                        }
                    });

                    for (auto &entry : bins) {
                        max_mean_cov = std::max(entry.second.mean_cov, max_mean_cov);
                    }

                    longest_path_len = std::max(longest_path_len, curr_len);
                }
            });

            max_mean_cov /= _bin_width;
        }

        std::unordered_set<pair<uint64_t, uint64_t>> edges_drawn;
        uint64_t gap_links_removed = 0;
        uint64_t total_links = 0;
        bool _color_path_names_background = args::get(color_path_names_background);

        graph.for_each_path_handle([&](const path_handle_t &path) {
            //std::cerr << "path " << as_integer(path) << " " << graph.get_path_name(path) << " " << path_r_f << " " << path_g_f << " " << path_b_f
            //          << " " << (int)path_r << " " << (int)path_g << " " << (int)path_b << std::endl;
            uint64_t path_rank = as_integer(path) - 1;
            if (path_layout_y[path_rank] >= 0){
                // use a sha256 to get a few bytes that we'll use for a color
                std::string path_name = graph.get_path_name(path);

    #ifdef debug_odgi_viz
                std::cerr << "path_name: " << path_name << std::endl;
    #endif

                bool is_aln = true;
                if (aln_mode) {
                    std::string::size_type n = path_name.find(aln_prefix);
                    if (n != 0) {
                        is_aln = false;
                    }
                }
                // use a sha256 to get a few bytes that we'll use for a color
                picosha2::byte_t hashed[picosha2::k_digest_size];
                if (_color_by_prefix) {
                    std::string path_name_prefix = prefix(path_name, path_name_prefix_separator);
                    picosha2::hash256(path_name_prefix.begin(), path_name_prefix.end(), hashed, hashed + picosha2::k_digest_size);
                } else {
                    picosha2::hash256(path_name.begin(), path_name.end(), hashed, hashed + picosha2::k_digest_size);
                }

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

                // Calculate the number or steps, the reverse steps and the length of the path if any of this information
                // is needed depending on the input arguments.
                uint64_t steps = 0;
                uint64_t rev = 0;
                uint64_t path_len_to_use = 0;
                std::map<uint64_t, algorithms::path_info_t> bins;
                if (is_aln) {
                    if (
                            _show_strands ||
                            (_change_darkness && !_longest_path) ||
                            (_binned_mode && (_color_by_mean_coverage || _color_by_mean_inversion_rate || _change_darkness))
                            ) {
                        handle_t h;
                        uint64_t hl, p;
                        bool is_rev;
                        graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                            h = graph.get_handle_of_step(occ);
                            is_rev = graph.get_is_reverse(h);
                            hl = graph.get_length(h);

                            if (_show_strands) {
                                ++steps;

                                rev += is_rev;
                            }

                            if (_change_darkness && !_longest_path) {
                                path_len_to_use += hl;
                            }

                            if (_binned_mode && (_color_by_mean_coverage || _color_by_mean_inversion_rate || _change_darkness)){
                                p = position_map[number_bool_packing::unpack_number(h)];
                                for (uint64_t k = 0; k < hl; ++k) {
                                    int64_t curr_bin = (p + k) / _bin_width + 1;

                                    ++bins[curr_bin].mean_cov;
                                    if (is_rev) {
                                        ++bins[curr_bin].mean_inv;
                                    }
                                }
                            }
                        });

                        if (_binned_mode && (_color_by_mean_coverage || _color_by_mean_inversion_rate || _change_darkness)) {
                            for (auto &entry : bins) {
                                auto &v = entry.second;
                                v.mean_inv /= (v.mean_cov ? v.mean_cov : 1);
                                v.mean_cov /= _bin_width;
                            }
                        }
                    }

                    if (_change_darkness && _longest_path){
                        path_len_to_use = longest_path_len;
                    }

                    if (_show_strands) {
                        float x = path_r_f;
                        path_r_f = (x + 0.5 * 9) / 10;
                        path_g_f = (x + 0.5 * 9) / 10;
                        path_b_f = (x + 0.5 * 9) / 10;
                        // check the path orientations
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
                    } else if (_change_darkness && _white_to_black) {
                        path_r = 220;
                        path_g = 220;
                        path_b = 220;
                    } else if (_binned_mode && _color_by_mean_coverage) {
                        path_r = 0;
                        path_g = 0;
                        path_b = 255;
                    } else if (_color_by_mean_inversion_rate) {
                        path_r = 255;
                        path_g = 0;
                        path_b = 0;
                    }
                }

                if (!(
                        is_aln && (( _change_darkness && _white_to_black) || _color_by_mean_inversion_rate || (_binned_mode && (_color_by_mean_coverage || _change_darkness)))
                        )) {
                    // brighten the color
                    float f = std::min(1.5, 1.0 / std::max(std::max(path_r_f, path_g_f), path_b_f));
                    path_r = (uint8_t) std::round(255 * std::min(path_r_f * f, (float) 1.0));
                    path_g = (uint8_t) std::round(255 * std::min(path_g_f * f, (float) 1.0));
                    path_b = (uint8_t) std::round(255 * std::min(path_b_f * f, (float) 1.0));
                }

                if (char_size >= 8){
                    uint8_t num_of_chars = min(path_name.length(), (uint64_t) max_num_of_chars);
                    bool path_name_too_long = path_name.length() > num_of_chars;

                    uint8_t ratio = char_size / 8;
                    uint8_t left_padding = max_num_of_chars - num_of_chars;

                    if (_color_path_names_background){
                        for (uint32_t x = left_padding * char_size; x <= max_num_of_chars * char_size; x++){
                            add_path_step(image_path_names, width_path_names, (float)(x + ratio) * (1 / scale_x), path_layout_y[path_rank], path_r, path_g, path_b);
                        }
                    }

                    uint64_t base_y = path_layout_y[path_rank] * pix_per_path + pix_per_path / 2 - char_size / 2;

                    for(uint16_t i = 0; i < num_of_chars; i++){
                        uint64_t base_x = (left_padding + i) * char_size;

                        auto cb = (i < num_of_chars - 1 || !path_name_too_long) ? font_5x8[path_name[i]] : font_5x8_special[TRAILING_DOTS];

                        write_character_in_matrix(
                                image_path_names, width_path_names, cb,
                                char_size,
                                base_x, base_y,
                                0, 0, 0
                        );
                    }
                }

                uint64_t curr_len = 0;
                double x = 1.0;
                if (_binned_mode) {
                    std::vector<std::pair<uint64_t, uint64_t>> links;
                    std::vector<uint64_t> bin_ids;
                    int64_t last_bin = 0; // flag meaning "null bin"

                    handle_t h;
                    uint64_t p, hl;

                    graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                        h = graph.get_handle_of_step(occ);
                        p = position_map[number_bool_packing::unpack_number(h)];
                        hl = graph.get_length(h);

                        // make contents for the bases in the node

                        uint64_t path_y = path_layout_y[path_rank];
                        for (uint64_t k = 0; k < hl; ++k) {
                            int64_t curr_bin = (p + k) / _bin_width + 1;

                            if (curr_bin != last_bin) {
                                bin_ids.push_back(curr_bin);

#ifdef debug_odgi_viz
                                std::cerr << "curr_bin: " << curr_bin << std::endl;
#endif

                                if (is_aln) {
                                    if (_change_darkness){
                                        uint64_t ii = bins[curr_bin].mean_inv > 0.5 ? (hl - k) : k;
                                        x = 1 - ( (float)(curr_len + ii) / (float)(path_len_to_use)) * 0.9;
                                    } else if (_color_by_mean_coverage) {
                                        x = bins[curr_bin].mean_cov / max_mean_cov;
                                    } else if (_color_by_mean_inversion_rate) {
                                        x = bins[curr_bin].mean_inv;
                                    }
                                }

                                if (curr_bin - 1 >= pangenomic_start_pos && curr_bin - 1 <= pangenomic_end_pos) {
                                    add_path_step(image, width, curr_bin - 1 - pangenomic_start_pos, path_y, (float)path_r * x, (float)path_g * x, (float)path_b * x);
                                }

                                if (std::abs(curr_bin - last_bin) > 1 || last_bin == 0) {
                                    // bin cross!
                                    links.push_back(std::make_pair(last_bin, curr_bin));
                                }
                            }

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
                                add_edge_from_positions(a, b, RGB_BIN_LINKS);
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
                            if (is_aln) {
                                if (_change_darkness){
                                    uint64_t ii = graph.get_is_reverse(h) ? (hl - i) : i;
                                    x = 1 - ((float)(curr_len + ii*scale_x) / (float)(path_len_to_use))*0.9;
                                } else if (_color_by_mean_inversion_rate) {
                                    if (graph.get_is_reverse(h)) {
                                        path_r = 255;
                                    } else {
                                        path_r = 0;
                                    }
                                };
                            }

                            if ((p + i) >= pangenomic_start_pos && (p + i) <= pangenomic_end_pos) {
                                add_path_step(image, width, p+i - pangenomic_start_pos, path_y, (float)path_r * x, (float)path_g * x, (float)path_b * x);
                            }
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
            }
        });

        if (args::get(drop_gap_links)) {
            std::cerr << std::setprecision(4) << "Gap links removed: " << (100.0 *  ((double)gap_links_removed / (double)total_links))
            << "%, that is " << gap_links_removed << " gap links (" << path_count << " path start links + "
            << path_count << " path end links + " << (gap_links_removed - path_count * 2) << " inner gap links) of "
            << total_links << " total links" << std::endl;
        }

        // trim horizontal and vertical spaces to fit
        uint64_t min_x = std::numeric_limits<uint64_t>::max();
        uint64_t max_x = std::numeric_limits<uint64_t>::min(); // 0
        uint64_t min_y = std::numeric_limits<uint64_t>::max();
        uint64_t max_y = std::numeric_limits<uint64_t>::min(); // 0
        for (uint64_t y = 0; y < height + path_space; ++y) {
            for (uint64_t x = 0; x < width; ++x) {
                uint8_t r = image[4 * width * y + 4 * x + 0];
                uint8_t g = image[4 * width * y + 4 * x + 1];
                uint8_t b = image[4 * width * y + 4 * x + 2];
                //uint8_t a = image[4 * width * y + 4 * x + 3];
                if (r != 255 || g != 255 || b != 255) {
                    min_x = std::min(min_x, x);
                    max_x = std::max(max_x, x);
                    min_y = std::min(min_y, y);
                    max_y = std::max(max_y, y);
                }
            }
        }

        // provide some default padding at the bottom, to clarify the edges
        max_y = std::min(path_space + height, max_y + bottom_padding);

        uint64_t crop_width = max_x - min_x + 1;
        uint64_t crop_height = max_y - min_y;

        if (char_size >= 8){
            crop_width += width_path_names;
        }

        /*std::cerr << "width " << width << std::endl;
        std::cerr << "height " << height << std::endl;
        std::cerr << "min and max x " << min_x << " " << max_x << std::endl;
        std::cerr << "min and max y " << min_y << " " << max_y << std::endl;
        std::cerr << "crop_width " << crop_width << std::endl;
        std::cerr << "crop_height " << crop_height << std::endl;*/

        std::vector<uint8_t> crop;
        crop.resize(crop_width * crop_height * 4, 255);
        for (uint64_t y = 0; y < crop_height; ++y) {
            for (uint64_t x = 0; x < crop_width; ++x) {
                for (uint8_t z = 0; z < 4; z++){
                    crop[4 * crop_width * y + 4 * x + z] = (char_size >= 8 && x < width_path_names) ?
                            image_path_names[4 * width_path_names * (y + min_y) + 4 * x + z] :
                            image[4 * width * (y + min_y) + 4 * (x - width_path_names + min_x) + z];
                }
            }
        }

        png::encodeOneStep(filename, crop, crop_width, crop_height);

        return 0;
    }

    static Subcommand odgi_viz("viz", "visualize the graph",PIPELINE, 3, main_viz);

}
