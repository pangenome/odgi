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
#include "utils.hpp"
#include "colorbrewer.hpp"
#include "split.hpp"

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
        return s.substr(0, s.find_last_of(c));
    }

    int main_viz(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi viz";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser("Visualize a variation graph in 1D.");
        args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
        args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
        args::ValueFlag<std::string> png_out_file(mandatory_opts, "FILE", "Write the visualization in PNG format to this *FILE*.", {'o', "out"});
        args::Group viz_opts(parser, "[ Visualization Options ]");
        args::ValueFlag<uint64_t> image_width(viz_opts, "N", "Set the width in pixels of the output image (default: 1500).", {'x', "width"});
        args::ValueFlag<uint64_t> image_height(viz_opts, "N", "Set the height in pixels of the output image (default: 500).", {'y', "height"});
        args::ValueFlag<uint64_t> path_height(viz_opts, "N", "The height in pixels for a path.", {'a', "path-height"});
        args::ValueFlag<uint64_t> path_x_pad(viz_opts, "N", "The padding in pixels on the x-axis for a path.", {'X', "path-x-padding"});
        args::Flag no_path_borders(viz_opts, "bool", "Don't show path borders.", {'n', "no-path-borders"});
        args::Flag black_path_borders(viz_opts, "bool", "Draw path borders in black (default is white).", {'b', "black-path-borders"});
        args::Flag pack_paths(viz_opts, "bool", "Pack all paths rather than displaying a single path per row.",{'R', "pack-paths"});
        args::ValueFlag<float> link_path_pieces(viz_opts, "FLOAT","Show thin links of this relative width to connect path pieces.",{'L', "link-path-pieces"});
        args::ValueFlag<std::string> alignment_prefix(viz_opts, "STRING","Apply alignment related visual motifs to paths which have this name"
                                                                         " prefix. It affects the [**-S, --show-strand**] and [**-d,"
                                                                         " –change-darkness**] options.",{'A', "alignment-prefix"});
        args::Flag show_strands(viz_opts, "bool","Use red and blue coloring to display forward and reverse alignments."
                                                 " This parameter can be set in combination with [**-A,"
                                                 " –alignment-prefix**=*STRING*].",{'S', "show-strand"});
        args::Flag color_by_mean_inversion_rate(viz_opts, "color-by-mean-inversion-rate", "Change the color respect to the node strandness (black for forward,"
                                                                                          " red for reverse); in binned mode (**-b, --binned-mode**), change the"
                                                                                          " color respect to the mean inversion rate of the path for each bin,"
                                                                                          " from black (no inversions) to red (bin mean inversion rate equals to"
                                                                                          " 1).", {'z', "color-by-mean-inversion-rate"});
        args::Flag color_by_uncalled_bases(viz_opts, "bool", "Change the color with respect to the uncalled bases of the path for each"
                                                             " bin, from black (no uncalled bases) to green (all uncalled bases).",
                                                              {'N', "color-by-uncalled-bases"});
        args::ValueFlag<char> color_by_prefix(viz_opts, "CHAR", "Color paths by their names looking at the prefix before the given"
                                                                 " character CHAR.",{'s', "color-by-prefix"});
        // TODO
        args::ValueFlag<std::string> _name_prefixes(viz_opts, "FILE", "Merge paths beginning with prefixes listed (one per line) in *FILE*.", {'M', "prefix-merges"});
        args::ValueFlag<std::string> _ignore_prefix(viz_opts, "PREFIX", "Ignore paths starting with the given *PREFIX*.", {'I', "ignore-prefix"});

        /// Range selection
        args::Group intervals_opts(parser, "[ Intervals Selection Options ]");
        args::ValueFlag<std::string> _nucleotide_range(intervals_opts, "STRING","Nucleotide range to visualize: ``STRING=[PATH:]start-end``. ``\\*-end``"
                                                                        " for ``[0,end]``; ``start-*`` for ``[start,pangenome_length]``. If no"
                                                                        " PATH is specified, the nucleotide positions refer to the pangenome’s"
                                                                        " sequence (i.e., the sequence obtained arranging all the graph’s node"
                                                                        " from left to right).",{'r', "path-range"});

        /// Paths selection
        args::Group path_selection_opts(parser, "[ Path Selection Options ]");
        args::ValueFlag<std::string> _path_names_file(path_selection_opts, "FILE", "List of paths to display in the specified order; the file must contain"
                                                                                   " one path name per line and a subset of all paths can be specified.", {'p', "paths-to-display"});

        /// Path names
        args::Group path_names_viz_opts(parser, "[ Path Names Viz Options ]");
        args::Flag hide_path_names(path_names_viz_opts, "bool", "Hide the path names on the left of the generated image.",{'H', "hide-path-names"});
        args::Flag color_path_names_background(path_names_viz_opts, "bool", "Color path names background with the same color as paths.",{'C', "color-path-names-background"});
        args::ValueFlag<size_t> _max_num_of_characters(path_names_viz_opts, "N", "Maximum number of characters to display for each path name (max 128"
                                                                                   " characters). The default value is *the length of the longest path"
                                                                                   " name* (up to 128 characters).",{'c', "max-num-of-characters"});

        /// Binned mode
        args::Group bin_opts(parser, "[ Binned Mode Options ]");
        /*
        args::Flag binned_mode(bin_opts, "binned-mode", "The variation graph is binned before its visualization. Each pixel in"
                                                        " the output image will correspond to a bin. For more information about"
                                                        " the binning process, please refer to odgi bin.", {'b', "binned-mode"});
        */
        args::ValueFlag<uint64_t> bin_width(bin_opts, "bp", "The bin width specifies the size of each bin in the binned mode. If it"
                                                            " is not specified, the bin width is calculated from the width in pixels"
                                                            " of the output image.r",{'w', "bin-width"});
        /*
        args::Flag drop_gap_links(bin_opts, "drop-gap-links", "We divide links into 2 classes:"
                                                              "1. The links which help to follow complex variations. They need to be"
                                                              " drawn, else one could not follow the sequence of a path."
                                                              "2. The links helping to follow simple variations. These links are called"
                                                              " **gap-links**. Such links solely connecting a path from left to right"
                                                              " may not be relevant to understand a path’s traversal through the"
                                                              " bins. Therefore, when this option is set, the gap-links are not drawn"
                                                              " in binned mode.", {'g', "no-gap-links"});
        */
        args::Flag color_by_mean_depth(bin_opts, "color-by-mean-depth", "Change the color with respect to the mean coverage of the path for each"
                                                                        " bin, using the colorbrewer palette specified in -B --colorbrewer-palette",
                                                                        {'m', "color-by-mean-depth"});
        args::ValueFlag<std::string> colorbrewer_palette(bin_opts, "SCHEME:N", "Use the colorbrewer palette specified by the given SCHEME, with the number of levels N."
                                                                               " Specifiy 'show' to see available palettes.",
                                                                               {'B', "colorbrewer-palette"});
        args::Flag no_grey_depth(bin_opts, "bool", "Use the colorbrewer palette for <0.5x and ~1x coverage bins."
                                 " By default, these bins are light and neutral grey.",
                                 {'G', "no-grey-depth"});

        /// Gradient mode
        args::Group grad_mode_opts(parser, "[ Gradient Mode Options ]");
        args::Flag change_darkness(grad_mode_opts, "change-darkness", "Change the color darkness based on nucleotide position in the path."
                                                                      " When it is used in binned mode, the mean inversion rate of the bin"
                                                                      " node is considered to set the color gradient starting position: when"
                                                                      " this rate is greater than 0.5, the bin is considered inverted, and the"
                                                                      " color gradient starts from the right-end of the bin. This parameter"
                                                                      " can be set in combination with [**-A,"
                                                                      " –alignment-prefix**=*STRING*].", {'d', "change-darkness"});
        args::Flag longest_path(grad_mode_opts, "longest-path", "Use the longest path length to change the color darkness.", {'l', "longest-path"});
        args::Flag white_to_black(grad_mode_opts, "white-to-black", "Change the color darkness from white (for the first nucleotide"
                                                                    " position) to black (for the last nucleotide position).", {'u', "white-to-black"});

		/// Compressed mode
		args::Group compressed_mode_opts(parser, "[ Compressed Mode Options ]");
		args::Flag compress(compressed_mode_opts, "compressed", "Compress the view vertically, summarizing the path coverage across all paths dis"
															  "playing the information using only one path 'COMPRESSED_MODE'. "
															  "A heatmap color-coding from https://colorbrewer2.org/#type=diverging&scheme=RdBu&n=11"
															  " is used. Alternatively, one can enter a colorbrewer palette via "
															  "-B, --colorbrewer-palette.", {'O', "compressed-mode"});

		args::Group threading(parser, "[ Threading ]");
		args::ValueFlag<uint64_t> nthreads(threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
		args::Group processing_info_opts(parser, "[ Processing Information ]");
		args::Flag _progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
        args::Group program_information(parser, "[ Program Information ]");
        args::HelpFlag help(program_information, "help", "Print a help message for odgi viz.", {'h', "help"});

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
                    << "[odgi::viz] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
                    << std::endl;
            return 1;
        }

        if (!png_out_file) {
            std::cerr
                    << "[odgi::viz] error: please specify an output file to where to store the PNG via -o=[FILE], --out=[FILE]."
                    << std::endl;
            return 1;
        }

        /*
        if (
                !args::get(binned_mode) &&
                ((args::get(bin_width) > 0) || args::get(drop_gap_links) ||
                args::get(color_by_mean_depth))
                ){
            std::cerr
                    << "[odgi::viz] error: please specify the -b/--binned-mode option to use the "
                       "-w/--bin_width, -g/--no-gap-links, and -m/--color-by-mean-depth "
                       "options."
                    << std::endl;
            return 1;
        }
        */

        if (!args::get(change_darkness) && (args::get(longest_path) || args::get(white_to_black))){
            std::cerr
                    << "[odgi::viz] error: please specify the -d/--change-darkness option to use the -l/--longest-path and -u/--white-to-black options."
                    << std::endl;
            return 1;
        }

        if ((args::get(color_by_prefix) != 0) + args::get(show_strands) + args::get(white_to_black) + args::get(color_by_mean_depth) + args::get(color_by_mean_inversion_rate) + args::get(color_by_uncalled_bases)  > 1) {
            std::cerr
                    << "[odgi::viz] error: please specify only one of the following options: "
                       "-s/--color-by-prefix, -S/--show-strand, -u/--white-to-black, "
                       "-m/--color-by-mean-depth, and -z/--color-by-mean-inversion."
                    << std::endl;
            return 1;
        }

        if (args::get(change_darkness) && (args::get(color_by_mean_depth) || args::get(color_by_mean_inversion_rate) || args::get(color_by_uncalled_bases))) {
            std::cerr
                    << "[odgi::viz] error: please specify the -d/--change-darkness option without specifying "
                       "-m/--color-by-mean-depth or -z/--color-by-mean-inversion or -N/--color-by-uncalled-bases."
                    << std::endl;
            return 1;
        }

        if (args::get(pack_paths) && !args::get(_path_names_file).empty()){
            std::cerr
                    << "[odgi::viz] error: please specify -R/--pack-paths or -p/--paths-to-display, not both."
                    << std::endl;
            return 1;
        }

        if (args::get(hide_path_names) && (args::get(color_path_names_background) || args::get(_max_num_of_characters))){
            std::cerr
                    << "[odgi::viz] error: please specify the -C/--color-path-names-background and -c/--max-num-of-characters "
                       "options without specifying -H/--hide-path-names." << std::endl;
            return 1;
        }

		const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

        //NOTE: this sample will overwrite the file or test.png without warning!
        //const char* filename = argc > 1 ? argv[1] : "test.png";
        if (args::get(png_out_file).empty()) {
            std::cerr << "[odgi::viz] error: output image required" << std::endl;
            return 1;
        }

        graph_t graph;
        assert(argc > 0);
        if (!args::get(dg_in_file).empty()) {
            const std::string infile = args::get(dg_in_file);
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
                utils::handle_gfa_odgi_input(infile, "viz", args::get(_progress), num_threads, graph);
            }
        }

        std::vector<uint64_t> position_map(graph.get_node_count() + 1);
        const uint64_t shift = number_bool_packing::unpack_number(graph.get_handle(graph.min_node_id()));
        uint64_t len = 0;
        {
            if (number_bool_packing::unpack_number(graph.get_handle(graph.max_node_id())) - shift >= graph.get_node_count()){
                std::cerr << "[odgi::viz] error: the node IDs are not compacted. Please run 'odgi sort' using -O, --optimize to optimize the graph." << std::endl;
                exit(1);
            }
            graph.for_each_handle([&](const handle_t &h) {
                position_map[number_bool_packing::unpack_number(h) - shift] = len;
                uint64_t hl = graph.get_length(h);
                len += hl;
#ifdef debug_odgi_viz
                std::cerr << "SEGMENT ID: " << graph.get_id(h) << " - " << as_integer(h) << " - index_in_position_map (" << (number_bool_packing::unpack_number(h) - shift) << ") = " << len << std::endl;
#endif
            });
            position_map[position_map.size() - 1] = len;
        }

        double pangenomic_start_pos = 0.0;
        double pangenomic_end_pos = (double) (len - 1);

        {
            std::string nucleotide_range = args::get(_nucleotide_range);
            if (!nucleotide_range.empty()) {
                const size_t foundLastColon = nucleotide_range.find_last_of(':');
                std::string path_name;
                if (foundLastColon != string::npos) {
                    path_name = nucleotide_range.substr(0, foundLastColon);

                    if (!graph.has_path(path_name)) {
                        std::cerr
                                << "[odgi::viz] error: please specify a valid path name."
                                << std::endl;
                        return 1;
                    }

                    nucleotide_range = nucleotide_range.substr(foundLastColon + 1);
                }

                const std::regex regex("-");
                const std::vector<std::string> splitted(
                        std::sregex_token_iterator(nucleotide_range.begin(), nucleotide_range.end(), regex, -1),
                        std::sregex_token_iterator()
                );

                if (splitted.size() != 2) {
                    std::cerr
                            << "[odgi::viz] error: please specify a valid nucleotide range: STRING=[PATH:]start-end."
                            << std::endl;
                    return 1;
                }

                if ((splitted[0] != "*" && !utils::is_number(splitted[0])) || (splitted[1] != "*" && !utils::is_number(splitted[1]))) {
                    std::cerr
                            << "[odgi::viz] error: please specify valid numbers for the nucleotide range."
                            << std::endl;
                    return 1;
                }

                auto get_path_length = [](const graph_t &graph, const path_handle_t &path_handle) {
                    uint64_t path_len = 0;
                    graph.for_each_step_in_path(path_handle, [&](const step_handle_t &s) {
                        path_len += graph.get_length(graph.get_handle_of_step(s));
                    });
                    return path_len;
                };
                const double path_len = (double) get_path_length(graph, graph.get_path_handle(path_name));

                if (splitted[0] == "*") {
                    pangenomic_start_pos = 0;
                } else {
                    pangenomic_start_pos = std::max(0.0, std::min(path_len - 1, stod(splitted[0])));
                }

                if (splitted[1] == "*") {
                    pangenomic_end_pos = (double)len - 1;
                } else {
                    pangenomic_end_pos = std::min(path_len, stod(splitted[1]));
                }

                //std::cerr << "input A: " << pangenomic_start_pos << std::endl;
                //std::cerr << "input B: " << pangenomic_end_pos << std::endl;

                if (!path_name.empty()) {
                    // Convert the nucleotide path range in a nucleotide pangenomic range
                    if(_progress) {
                        std::cerr
                        << "[odgi::viz] Path range to pangenomic range conversion."
                        << std::endl;
                    }

                    double new_pangenomic_start_pos = (double) (len - 1);
                    double new_pangenomic_end_pos = 0;

                    const path_handle_t path_handle = graph.get_path_handle(path_name);

                    uint64_t nt_position_in_path = 0;

                    graph.for_each_step_in_path(path_handle, [&](const step_handle_t &occ) {
                        const handle_t h = graph.get_handle_of_step(occ);
                        uint64_t h_pan_pos = position_map[number_bool_packing::unpack_number(h) - shift];
                        uint64_t h_len = graph.get_length(h);

                        //Todo dumb implementation: improve with the math later
                        do {
                            if (nt_position_in_path >= pangenomic_start_pos && nt_position_in_path <= pangenomic_end_pos) {
                                new_pangenomic_start_pos = std::min(new_pangenomic_start_pos, (double) h_pan_pos);
                                new_pangenomic_end_pos = std::max(new_pangenomic_end_pos, (double) h_pan_pos);
                            }

                            ++nt_position_in_path;
                            ++h_pan_pos;
                        } while(--h_len != 0);
                    });

                    //std::cerr << "A: " << new_pangenomic_start_pos << std::endl;
                    //std::cerr << "B: " << new_pangenomic_end_pos << std::endl;

                    pangenomic_start_pos = new_pangenomic_start_pos;
                    pangenomic_end_pos = new_pangenomic_end_pos;
                }

                if (pangenomic_start_pos >= pangenomic_end_pos) {
                    std::cerr
                            << "[odgi::viz] error: please specify a start position less than the end position."
                            << std::endl;
                    return 1;
                }

                if(_progress) {
                    std::cerr
                    << "[odgi::viz] Visualizing the graph in the pangenomic range [" << pangenomic_start_pos << ", " << pangenomic_end_pos << "]"
                    << std::endl;
                }
            }
        }

        const size_t max_num_of_characters = args::get(_max_num_of_characters) > 1 ? min(args::get(_max_num_of_characters), (size_t) PATH_NAMES_MAX_NUM_OF_CHARACTERS) : PATH_NAMES_MAX_NUM_OF_CHARACTERS;

        uint64_t path_count = graph.get_path_count();
        const uint64_t pix_per_path = args::get(path_height) ? args::get(path_height) : 10;
        const uint64_t pix_per_link = std::max((uint64_t) 1, (uint64_t) std::round(args::get(link_path_pieces) * pix_per_path));
        const uint64_t link_pix_y = pix_per_path / 2 - pix_per_link / 2;
        const uint64_t path_padding = args::get(path_x_pad);
        const uint64_t bottom_padding = 5;

        const uint64_t len_to_visualize = pangenomic_end_pos - pangenomic_start_pos + 1;

        // the math here only works if the image size is greater than or equal to the graph length
        uint64_t width = std::min(len_to_visualize, (args::get(image_width) ? args::get(image_width) : 1500));
        const uint64_t height = std::min(len_to_visualize, (args::get(image_height) ? args::get(image_height) : 500) + bottom_padding);

        double scale_x = (double) width / (double) len_to_visualize;
        const double scale_y = (double) height / (double) len_to_visualize;

        double _bin_width = args::get(bin_width);
        bool _binned_mode = true; //args::get(binned_mode);
        if (_binned_mode){
            if (_bin_width == 0){
                _bin_width = 1.0 / scale_x; // It can be a double value.
            }else{
                width = len_to_visualize / _bin_width;// + (len_to_visualize % bin_width ? 1 : 0);
            }

            pangenomic_start_pos /= _bin_width;
            pangenomic_end_pos /= _bin_width;

            // Each pixel corresponds to a bin
            scale_x = 1; //scale_x*bin_width;

            if(_progress) {
                std::cerr << "[odgi::viz] Binned mode" << std::endl;
                std::cerr << "[odgi::viz] bin width: " << _bin_width << std::endl;
                std::cerr << "[odgi::viz] image width: " << width << std::endl;
            }
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

        std::vector<int64_t> path_group;
        std::vector<std::string> prefixes;
        bool group_paths = false;
        if (_name_prefixes) {
            std::vector<std::string> prefixes_tmp;
            group_paths = true;
            auto& name_prefixes = args::get(_name_prefixes);
            ska::flat_hash_set<std::string> group_names;
            if (!name_prefixes.empty()){
                std::ifstream prefixes_in(name_prefixes);
                std::string line;
                while (std::getline(prefixes_in, line)) {
                    if (!line.empty()){
                        if (group_names.count(line)) {
                            std::cerr << "[odgi::viz] duplicate group name found " << line << std::endl;
                        } else {
                            prefixes_tmp.push_back(line);
                            group_names.insert(line);
                        }
                    }
                }
            }

            path_group.resize(path_count, -1);
            graph.for_each_path_handle(
                [&](const path_handle_t &path) {
                    uint64_t path_rank = as_integer(path) - 1;
                    // XXX quadratic
                    // could be O(N log N) if we use binary search on a sorted list

                    auto path_name = graph.get_path_name(path);
                    uint64_t group_idx = 0;
                    bool found = false;

                    // Search first in already validated prefixes
                    for (auto& prefix : prefixes) {
                        // check if it matches the start
                        if (path_name.find(prefix) == 0) {
                            found = true;
                            path_group[path_rank] = group_idx;
                            break;
                        }
                        ++group_idx;
                    }

                    // If not found in the validated prefixes...
                    if (!found) {
                        // ... search in all read prefixes
                        for (auto& prefix : prefixes_tmp) {
                            // check if it matches the start
                            if (path_name.find(prefix) == 0) {
                                prefixes.push_back(prefix); // Add into the validated prefixes
                                path_group[path_rank] = group_idx;
                                break;
                            }
                        }
                    }
                });

            if (_progress) {
                std::cerr << "[odgi::viz] Read " << prefixes.size() << " valid prefixes." << std::endl;
            }
        }

        const std::string ignore_prefix = _ignore_prefix ? args::get(_ignore_prefix) : "";

        const std::string compressed_path_name = "COMPRESSED_MODE";

        auto get_path_display_name =
            [&](const path_handle_t& p) {
					if (!ignore_prefix.empty()
						&& graph.get_path_name(p).find(ignore_prefix) == 0) {
						return std::string();
					} else {
						if (group_paths) {
							return prefixes[path_group[as_integer(p) - 1]];
						} else {
							return graph.get_path_name(p);
						}
					}
            };

        auto get_path_idx =
            [&](const path_handle_t& p) {
			if (compress) {
				return (int64_t) 0;
			} else {
				if (!ignore_prefix.empty()
					&& graph.get_path_name(p).find(ignore_prefix) == 0) {
					return (int64_t) - 1;
				} else {
					if (group_paths) {
						return (int64_t) path_group[as_integer(p) - 1];
					} else {
						return (int64_t)(as_integer(p) - 1);
					}
				}
			}
            };

        // map from path id to its starting y position
        //hash_map<uint64_t, uint64_t> path_layout_y;
        std::vector<int64_t> path_layout_y;
        path_layout_y.resize(path_count, -1);
        if (!args::get(pack_paths)) {
            std::string path_names = args::get(_path_names_file);
            if (!path_names.empty()){
                std::ifstream path_names_in(path_names);

                uint64_t rank_for_visualization = 0;
                uint64_t num_of_paths_in_file = 0;

                std::string line;
                while (std::getline(path_names_in, line)) {
                    if (!line.empty()){
                        if (graph.has_path(line)){
                            // todo here we need to do our grouping
                            int64_t path_rank = get_path_idx(graph.get_path_handle(line));
                            if (path_rank >= 0 && path_layout_y[path_rank] < 0) {
                                path_layout_y[path_rank] = rank_for_visualization++;
                                //rank_for_visualization++;
                            } else if (!group_paths) {
                                std::cerr << "[odgi::viz] error: in the path list there are duplicated path names." << std::endl;
                                exit(1);
                            }
                        }

                        num_of_paths_in_file++;
                    }
                }

                path_names_in.close();

                if (_progress){
                    std::cerr << "[odgi::viz] Found " << rank_for_visualization << "/" << num_of_paths_in_file << " paths to display." << std::endl;
                }

                if (rank_for_visualization == 0){
                    std::cerr << "[odgi::viz] error: no path to display." << std::endl;
                    exit(1);
                }

                path_count = rank_for_visualization;
            }else{
                uint64_t rank_for_visualization = 0;
                for (uint64_t i = 0; i < path_count; ++i) {
                    // todo here we need to do our grouping
                    int64_t path_rank = get_path_idx(as_path_handle(i+1));
                    //std::cerr << path_rank << " - " << path_layout_y[path_rank] << " - " << graph.get_path_name(as_path_handle(i+1)) << " --> " << rank_for_visualization << std::endl;
                    if (path_rank >= 0 && path_layout_y[path_rank] < 0) {
                        path_layout_y[path_rank] = rank_for_visualization++;
                    }
                }
                path_count = rank_for_visualization;
            }
        } else { // pack the paths
            if (group_paths) {
                std::cerr << "[odgi::viz] error: not possible to combine path packing and grouping" << std::endl;
                exit(1);
            }
            // buffer to record layout bounds
            std::vector<bool> path_layout_buf;
            path_layout_buf.resize(path_count * width);
            std::vector<path_handle_t> path_order = algorithms::id_ordered_paths(graph);
            for (auto &path : path_order) {
                // get the block which this path covers
                uint64_t min_x = len_to_visualize;
                uint64_t max_x = std::numeric_limits<uint64_t>::min(); // 0
                graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                    handle_t h = graph.get_handle_of_step(occ);
                    uint64_t p = position_map[number_bool_packing::unpack_number(h) - shift];

                    if (p >= pangenomic_start_pos && p <= pangenomic_end_pos) {
                        min_x = std::min(min_x, (uint64_t) (p - pangenomic_start_pos));
                        max_x = std::max(max_x, (uint64_t) (p + graph.get_length(h) - pangenomic_start_pos));
                    }
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

        const uint64_t path_space = path_count * pix_per_path;

        std::vector<uint8_t> image;
        image.resize(width * (height + path_space) * 4, 255);

        std::vector<uint8_t> image_path_names;
        if (!args::get(hide_path_names) && !args::get(pack_paths) && pix_per_path >= 8) {
            size_t _max_num_of_chars = std::numeric_limits<size_t>::min();

            if (group_paths) {
                for (auto& prefix : prefixes) {
                    _max_num_of_chars = max(_max_num_of_chars, prefix.length());
                }
            } else {
                graph.for_each_path_handle(
                        [&](const path_handle_t &path) {
                            int64_t path_rank = get_path_idx(path);
                            if (path_rank >= 0 && path_layout_y[path_rank] >= 0){
                                _max_num_of_chars = max(_max_num_of_chars, graph.get_path_name(path).length());
                            }
                        });
            }

            // In this way, compress_path_name will take up the same amount of space as the paths in uncompressed mode.
            // This makes it easier to compare/view the two viewing modes, easing interpretation.
            if (compress) {
                _max_num_of_chars = max(_max_num_of_chars, compressed_path_name.length());
            }

            max_num_of_chars = min(_max_num_of_chars, max_num_of_characters);

            char_size = min((uint16_t)((pix_per_path / 8) * 8), (uint16_t) PATH_NAMES_MAX_CHARACTER_SIZE);

            width_path_names = max_num_of_chars * char_size + char_size / 2;
            image_path_names.resize(width_path_names * (height + path_space) * 4, 255);
        }

        if (width_path_names + width > 50000){
            std::cerr
                    << "[odgi::viz] warning: you are going to create a big image (width > 50000 pixels)."
                    << std::endl;
        }

        const bool aln_mode = !args::get(alignment_prefix).empty();
        std::string aln_prefix;
        if (aln_mode) {
            aln_prefix = args::get(alignment_prefix);
        }

        char path_name_prefix_separator = '\0';
        if (color_by_prefix) {
            path_name_prefix_separator = args::get(color_by_prefix);
        }

        auto add_point = [&](const double &_x, const double &_y,
                             const uint8_t &_r, const uint8_t &_g, const uint8_t &_b) {
            uint64_t x = std::min((uint64_t) std::round(_x * scale_x), width - 1);
            uint64_t y = std::min((uint64_t) std::round(_y * scale_y), height - 1) + path_space;
            image[4 * width * y + 4 * x + 0] = _r;
            image[4 * width * y + 4 * x + 1] = _g;
            image[4 * width * y + 4 * x + 2] = _b;
            image[4 * width * y + 4 * x + 3] = 255;
        };

        auto add_edge_from_positions = [&](double a, const double b, uint8_t rgb) {
#ifdef debug_odgi_viz
            std::cerr << "Edge displayed" << std::endl;
            std::cerr << a << " --> " << b << std::endl;
#endif
            // Show the link if at least one of its 2 extremes is visible
            if (a >= pangenomic_start_pos && a <= pangenomic_end_pos || b >= pangenomic_start_pos && b <= pangenomic_end_pos) {
                // In binned mode, the Links have to be tall to be visible; in standard mode, _bin_width is 1, so nothing changes here
                const uint64_t dist = (b - a) * _bin_width;
                
                double i = 0.0;

                for (; i < dist; i += 1.0 / scale_y) {
                    if (a >= pangenomic_start_pos && a <= pangenomic_end_pos) {
                        add_point(a - pangenomic_start_pos, i, rgb, rgb, rgb);
                    }
                }

                while (a <= b) {
                    if (a >= pangenomic_start_pos && a <= pangenomic_end_pos) {
                        add_point(a - pangenomic_start_pos, i, rgb, rgb, rgb);
                    }
                    a += 1.0 / scale_x;
                }
                if (b >= pangenomic_start_pos && b <= pangenomic_end_pos) {
                    for (double j = 0.0; j < dist; j += 1.0 / scale_y) {
                        add_point(b - pangenomic_start_pos, j, rgb, rgb, rgb);
                    }
                }
            }
        };

        auto add_edge_from_handles = [&](const handle_t& h, const handle_t& o) {
            // map into our bins
            const uint64_t index_h = number_bool_packing::unpack_number(h) + !number_bool_packing::unpack_bit(h) - shift;
            const uint64_t index_o = number_bool_packing::unpack_number(o) + number_bool_packing::unpack_bit(o) - shift;
            const uint64_t h_pos = position_map[index_h] / _bin_width;
            const uint64_t o_pos = position_map[index_o] / _bin_width;

            // The last node has to be treated differently, as it has no following node, and its outgoing links would start outside the image
            const double x_shift = (index_h == position_map.size() - 1) || (index_o == position_map.size() - 1) ? (1.0 / _bin_width) : 0.0;

            const uint64_t a = std::min(h_pos, o_pos);
            const uint64_t b = std::max(h_pos, o_pos) >= x_shift ? (std::max(h_pos, o_pos) - x_shift) : 0;

#ifdef debug_odgi_viz
            std::cerr << graph.get_id(h) << " (" << number_bool_packing::unpack_bit(h) << ") --> " << graph.get_id(o) << " (" << number_bool_packing::unpack_bit(o) << ") " << std::endl;
            std::cerr << "edge " << a << " --> " << b << std::endl;
#endif

            add_edge_from_positions(std::min(a, b), std::max(a, b), 0);
        };

        {
            std::function<bool(const handle_t)> is_a_handle_to_hide;
			if (compress) {
				is_a_handle_to_hide = [&](const handle_t &h) {
					return false;
				};
			} else {
				if (path_count < graph.get_path_count()) {
					is_a_handle_to_hide = [&](const handle_t &h) {
						return graph.for_each_step_on_handle(h, [&](const step_handle_t &step) {
							auto path_handle = graph.get_path_handle_of_step(step);
							if (path_layout_y[as_integer(path_handle) - 1] >= 0) {
								return false;
							}

							return true;
						});
					};
				} else {
					is_a_handle_to_hide = [&](const handle_t &h) {
						return false;
					};
				}
			}

            /* FIXME Can we remove this?
            if (_binned_mode){
                graph.for_each_handle([&](const handle_t &h) {
                    if (!is_a_handle_to_hide(h)){
                        uint64_t p = position_map[number_bool_packing::unpack_number(h) - shift];
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
            }
            */
            graph.for_each_handle([&](const handle_t &h) {
                if (!is_a_handle_to_hide(h)){
                    uint64_t p = position_map[number_bool_packing::unpack_number(h) - shift];
                    uint64_t hl = graph.get_length(h);
                    // make contents for the bases in the node
                    for (double i = 0.0; i < hl; i += 1.0 / scale_x) {
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

        bool _no_path_borders = args::get(no_path_borders);
        bool _black_border = args::get(black_path_borders);
        auto add_path_step = [&](std::vector<uint8_t> &img, uint64_t width_img,
                const double &_x, const double &_y,
                const uint8_t &_r, const uint8_t &_g, const uint8_t &_b) {
            const uint64_t x = std::min((uint64_t) std::round(_x * scale_x), width - 1);
            const uint64_t t = _y * pix_per_path;
            uint64_t y = t;
            if (_no_path_borders || pix_per_path < 3) {
                const uint64_t s = t + pix_per_path;
                for ( ; y < s; ++y) {
                    img[4 * width_img * y + 4 * x + 0] = _r;
                    img[4 * width_img * y + 4 * x + 1] = _g;
                    img[4 * width_img * y + 4 * x + 2] = _b;
                    img[4 * width_img * y + 4 * x + 3] = 255;
                }
            } else {
                const uint64_t s = t + pix_per_path - 1;
                for ( ; y < s; ++y) {
                    img[4 * width_img * y + 4 * x + 0] = _r;
                    img[4 * width_img * y + 4 * x + 1] = _g;
                    img[4 * width_img * y + 4 * x + 2] = _b;
                    img[4 * width_img * y + 4 * x + 3] = 255;
                }
                if (_black_border) {
                    for ( ; y < s + 1; ++y) {
                        img[4 * width_img * y + 4 * x + 0] = 0;
                        img[4 * width_img * y + 4 * x + 1] = 0;
                        img[4 * width_img * y + 4 * x + 2] = 0;
                        img[4 * width_img * y + 4 * x + 3] = 255;
                    }
                }
            }
        };

        auto add_path_link = [&](const double &_x, const double &_y,
                                 const uint8_t &_r, const uint8_t &_g, const uint8_t &_b) {
            const uint64_t x = std::min((uint64_t) std::round(_x * scale_x), width - 1);
            const uint64_t t = _y * pix_per_path + link_pix_y;
            const uint64_t s = t + pix_per_link;
            for (uint64_t y = t; y < s; ++y) {
                image[4 * width * y + 4 * x + 0] = _r;
                image[4 * width * y + 4 * x + 1] = _g;
                image[4 * width * y + 4 * x + 2] = _b;
                image[4 * width * y + 4 * x + 3] = 255;
            }
        };

        const bool _show_strands = args::get(show_strands);

        const bool _change_darkness = args::get(change_darkness);
        const bool _longest_path = args::get(longest_path);
        const bool _white_to_black = args::get(white_to_black);

        const bool _color_by_mean_depth = args::get(color_by_mean_depth);
        const bool _color_by_mean_inversion_rate = args::get(color_by_mean_inversion_rate);
        const bool _color_by_uncalled_bases = args::get(color_by_uncalled_bases);

        uint64_t longest_path_len = 0;
        if ((_change_darkness && _longest_path) || (_binned_mode && _color_by_mean_depth)){
            graph.for_each_path_handle([&](const path_handle_t &path) {
                int64_t path_rank = get_path_idx(path);
                if (path_rank >= 0 && path_layout_y[path_rank] >= 0){
                    uint64_t curr_len = 0, p, hl;
                    handle_t h;
                    graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                        h = graph.get_handle_of_step(occ);
                        hl = graph.get_length(h);

                        curr_len += hl;
                    });

                    longest_path_len = std::max(longest_path_len, curr_len);
                }
            });
        }

        uint64_t gap_links_removed = 0;
        uint64_t total_links = 0;
        const bool _color_path_names_background = args::get(color_path_names_background);

		// Compressed-Mode part starts here :)
		if (compress) {
			std::map <uint64_t, algorithms::path_info_t> bins;
			graph.for_each_path_handle([&](const path_handle_t &path) {
				graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
					handle_t h = graph.get_handle_of_step(occ);
					uint64_t hl = graph.get_length(h);

					uint64_t p = position_map[number_bool_packing::unpack_number(h) - shift];
					for (uint64_t k = 0; k < hl; ++k) {
						int64_t curr_bin = (p + k) / _bin_width + 1;
						++bins[curr_bin].mean_depth;
					}
				});
			});

			/// path name part

			uint8_t path_r = 255;
			uint8_t path_g = 255;
			uint8_t path_b = 255;

			uint64_t path_rank = 0;

            const uint8_t num_of_chars = min(compressed_path_name.length(), (size_t)max_num_of_chars);
            const bool path_name_too_long = compressed_path_name.length() > num_of_chars;

            const uint8_t left_padding = max_num_of_chars - num_of_chars;

			// TODO Do we want this functionality?
			// uint8_t ratio = char_size / 8;
			/*
			if (_color_path_names_background) {
				for (uint32_t x = left_padding * char_size; x <= max_num_of_chars * char_size; x++) {
					add_path_step(image_path_names, width_path_names,
								  (double) (x + ratio) * (1.0 / scale_x), path_layout_y[path_rank], path_r,
								  path_g, path_b);
				}
			}
			 */
            const uint64_t base_y = path_layout_y[path_rank] * pix_per_path + pix_per_path / 2 - char_size / 2;

			for (uint16_t i = 0; i < num_of_chars; i++) {
                const uint64_t base_x = (left_padding + i) * char_size;

				auto cb = (i < num_of_chars - 1 || !path_name_too_long) ? font_5x8[compressed_path_name[i]]
																		: font_5x8_special[TRAILING_DOTS];

				write_character_in_matrix(
						image_path_names, width_path_names, cb,
						char_size,
						base_x, base_y,
						0, 0, 0
				);
			}
			/// end path name part

			double x = 1.0;
			colorbrewer::palette_t cov_colors;
			std::vector<double> cov_cuts;
			// Let the user enter the color palette
			if (colorbrewer_palette) {
				const auto parts = split(args::get(colorbrewer_palette), ':');
				cov_colors = colorbrewer::get_palette(parts.front(), std::stoi(parts.back()));
			} else {
				// we also have a default color palette https://colorbrewer2.org/#type=diverging&scheme=RdBu&n=11
				cov_colors = colorbrewer::get_palette("RdBu", 11);
			}
			uint64_t i = 0;
			cov_cuts.resize(cov_colors.size());
			double depth = 0.5;
			for (auto &color: cov_colors) {
				cov_cuts[i++] = depth;
				depth += 1;
			}

			for (auto &entry: bins) {
				auto &sec = entry.second;
				auto &curr_bin = entry.first;
				sec.mean_depth /= _bin_width;
				// std::cerr << "MEAN DEPTH OF BIN: " << v.mean_depth << std::endl;
				auto &mean_depth = sec.mean_depth;
				uint64_t j = 0;
				for (; j < cov_cuts.size(); ++j) {
					if (mean_depth <= cov_cuts[j]) {
						auto &v = cov_colors[j];
						path_r = v.red;
						path_g = v.green;
						path_b = v.blue;
						break;
					}
				}
				// take the max color
				if (j == cov_cuts.size()) {
					auto &v = cov_colors[j - 1];
					path_r = v.red;
					path_g = v.green;
					path_b = v.blue;
				}
				uint64_t path_y = path_layout_y[path_rank];
				add_path_step(image, width, curr_bin - 1 - pangenomic_start_pos, path_y,
							  (float) path_r * x, (float) path_g * x, (float) path_b * x);
			}
			/// end compressed-mode

			/// default case:
		} else {

			graph.for_each_path_handle([&](const path_handle_t &path) {
				int64_t path_rank = get_path_idx(path);
				//std::cerr << graph.get_path_name(path) << " -> " << path_rank << std::endl;
				if (path_rank >= 0 && path_layout_y[path_rank] >= 0) {
					// use a sha256 to get a few bytes that we'll use for a color
					std::string path_name = get_path_display_name(path);

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
					if (color_by_prefix) {
						std::string path_name_prefix = prefix(path_name, path_name_prefix_separator);
						picosha2::hash256(path_name_prefix.begin(), path_name_prefix.end(), hashed,
										  hashed + picosha2::k_digest_size);
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
					std::map <uint64_t, algorithms::path_info_t> bins;
					if (is_aln) {
						if (
								_show_strands ||
								(_change_darkness && !_longest_path) ||
								(_binned_mode &&
								 (_color_by_mean_depth || _color_by_mean_inversion_rate || _change_darkness ||
								  _color_by_uncalled_bases))
								) {
							handle_t h;
							uint64_t hl, p;
							bool is_rev;
							uint64_t num_uncalled_bases;
							graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
								h = graph.get_handle_of_step(occ);
								is_rev = graph.get_is_reverse(h);
								hl = graph.get_length(h);

								if (_color_by_uncalled_bases) {
									num_uncalled_bases = 0;
									for (auto c: graph.get_sequence(h)) {
										if (c == 'N' || c == 'n') {
											num_uncalled_bases++;
										}
									}
								}

								if (_show_strands) {
									++steps;

									rev += is_rev;
								}

								if (_change_darkness && !_longest_path) {
									path_len_to_use += hl;
								}

								if (_binned_mode &&
									(_color_by_mean_depth || _color_by_mean_inversion_rate || _change_darkness)) {
									p = position_map[number_bool_packing::unpack_number(h) - shift];
									for (uint64_t k = 0; k < hl; ++k) {
										int64_t curr_bin = (p + k) / _bin_width + 1;

										++bins[curr_bin].mean_depth;
										if (is_rev) {
											++bins[curr_bin].mean_inv;
										}
									}
								} else if (_binned_mode && _color_by_uncalled_bases) {
									p = position_map[number_bool_packing::unpack_number(h) - shift];
									for (uint64_t k = 0; k < hl; ++k) {
										int64_t curr_bin = (p + k) / _bin_width + 1;

										// Use the `mean_depth` field as 'mean_Ns`
										bins[curr_bin].mean_depth += num_uncalled_bases;
									}
								}
							});

							if (_binned_mode &&
								(_color_by_mean_depth || _color_by_mean_inversion_rate || _change_darkness)) {
								for (auto &entry: bins) {
									auto &v = entry.second;
									v.mean_inv /= (v.mean_depth ? v.mean_depth : 1);
									v.mean_depth /= _bin_width;
								}
							} else if (_binned_mode && _color_by_uncalled_bases) {
								for (auto &entry: bins) {
									auto &v = entry.second;
									v.mean_depth /= _bin_width;
								}
							}
						}

						if (_change_darkness && _longest_path) {
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
						} else if (_color_by_mean_inversion_rate) {
							path_r = 255;
							path_g = 0;
							path_b = 0;
						} else if (_color_by_uncalled_bases) {
							path_r = 0;
							path_g = 255;
							path_b = 0;
						}
					}

					if (!(
							is_aln && ((_change_darkness && _white_to_black) || _color_by_mean_inversion_rate ||
									   (_binned_mode &&
										(_color_by_mean_depth || _change_darkness || _color_by_uncalled_bases)))
					)) {
						// brighten the color
						float f = std::min(1.5, 1.0 / std::max(std::max(path_r_f, path_g_f), path_b_f));
						path_r = (uint8_t) std::round(255 * std::min(path_r_f * f, (float) 1.0));
						path_g = (uint8_t) std::round(255 * std::min(path_g_f * f, (float) 1.0));
						path_b = (uint8_t) std::round(255 * std::min(path_b_f * f, (float) 1.0));
					}

					if (char_size >= 8) {
                        const uint8_t num_of_chars = min(path_name.length(), (size_t) max_num_of_chars);
                        const bool path_name_too_long = path_name.length() > num_of_chars;

                        const uint8_t ratio = char_size / 8;
						const uint8_t left_padding = max_num_of_chars - num_of_chars;

						if (_color_path_names_background) {
							for (uint32_t x = left_padding * char_size; x <= max_num_of_chars * char_size; x++) {
								add_path_step(image_path_names, width_path_names,
											  (double) (x + ratio) * (1.0 / scale_x), path_layout_y[path_rank], path_r,
											  path_g, path_b);
							}
						}

                        const uint64_t base_y = path_layout_y[path_rank] * pix_per_path + pix_per_path / 2 - char_size / 2;

						for (uint16_t i = 0; i < num_of_chars; i++) {
							uint64_t base_x = (left_padding + i) * char_size;

							auto cb = (i < num_of_chars - 1 || !path_name_too_long) ? font_5x8[path_name[i]]
																					: font_5x8_special[TRAILING_DOTS];

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
						colorbrewer::palette_t cov_colors;
						std::vector<double> cov_cuts;
						if (_color_by_mean_depth) {
							if (colorbrewer_palette) {
								const auto parts = split(args::get(colorbrewer_palette), ':');
								cov_colors = colorbrewer::get_palette(parts.front(), std::stoi(parts.back()));
							} else {
								cov_colors = colorbrewer::get_palette("Spectral", 11);
							}
							if (!args::get(no_grey_depth)) {
								std::reverse(cov_colors.begin(),
											 cov_colors.end());
								cov_colors.push_back({128, 128, 128});
								cov_colors.push_back({196, 196, 196});
								std::reverse(cov_colors.begin(),
											 cov_colors.end());
							}
							uint64_t i = 0;
							cov_cuts.resize(cov_colors.size());
							double depth = 0.5;
							for (auto &color: cov_colors) {
								cov_cuts[i++] = depth;
								depth += 1;
							}
						}

						std::vector <std::pair<uint64_t, uint64_t>> links;
						std::vector <uint64_t> bin_ids;
						int64_t last_bin = 0; // flag meaning "null bin"

						handle_t h;
						uint64_t p, hl;

						graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
							h = graph.get_handle_of_step(occ);
							p = position_map[number_bool_packing::unpack_number(h) - shift];
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
										if (_change_darkness) {
											uint64_t ii = bins[curr_bin].mean_inv > 0.5 ? (hl - k) : k;
											x = 1.0 - ((double) (curr_len + ii) / (double) (path_len_to_use)) * 0.9;
										} else if (_color_by_mean_depth) {
											auto &mean_depth = bins[curr_bin].mean_depth;
											uint64_t j = 0;
											for (; j < cov_cuts.size(); ++j) {
												if (mean_depth <= cov_cuts[j]) {
													auto &v = cov_colors[j];
													path_r = v.red;
													path_g = v.green;
													path_b = v.blue;
													break;
												}
											}
											// take the max color
											if (j == cov_cuts.size()) {
												auto &v = cov_colors[j - 1];
												path_r = v.red;
												path_g = v.green;
												path_b = v.blue;
											}
										} else if (_color_by_mean_inversion_rate) {
											x = bins[curr_bin].mean_inv;
										} else if (_color_by_uncalled_bases) {
											x = bins[curr_bin].mean_depth;
										}
									}

									if (curr_bin - 1 >= pangenomic_start_pos && curr_bin - 1 <= pangenomic_end_pos) {
										add_path_step(image, width, curr_bin - 1 - pangenomic_start_pos, path_y,
													  (float) path_r * x, (float) path_g * x, (float) path_b * x);
									}

								}

								last_bin = curr_bin;
							}

							curr_len += hl;
						});

					} else {
						/// Loop over all the steps along a path, from first through last and draw them
						graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
							handle_t h = graph.get_handle_of_step(occ);
							uint64_t p = position_map[number_bool_packing::unpack_number(h) - shift];
							uint64_t hl = graph.get_length(h);
							// make contects for the bases in the node
							uint64_t path_y = path_layout_y[path_rank];
							for (uint64_t i = 0; i < hl; i += 1 / scale_x) {
								if (is_aln) {
									if (_change_darkness) {
										uint64_t ii = graph.get_is_reverse(h) ? (hl - i) : i;
										x = 1.0 -
											((double) (curr_len + ii * scale_x) / (double) (path_len_to_use)) * 0.9;
									} else if (_color_by_mean_inversion_rate) {
										if (graph.get_is_reverse(h)) {
											path_r = 255;
										} else {
											path_r = 0;
										}
									};
								}

								if ((p + i) >= pangenomic_start_pos && (p + i) <= pangenomic_end_pos) {
									add_path_step(image, width, p + i - pangenomic_start_pos, path_y,
												  (float) path_r * x, (float) path_g * x, (float) path_b * x);
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
							uint64_t p = position_map[number_bool_packing::unpack_number(h) - shift];
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
				//add_point(curr_bin - 1 - pangenomic_start_pos, 0, RGB_BIN_LINKS, RGB_BIN_LINKS, RGB_BIN_LINKS);
			});

		}

        /*
        if (args::get(drop_gap_links)) {
            std::cerr << std::setprecision(4) << "[odgi::viz] Gap links removed: " << (100.0 *  ((double)gap_links_removed / (double)total_links))
            << "%, that is " << gap_links_removed << " gap links (" << path_count << " path start links + "
            << path_count << " path end links + " << (gap_links_removed - path_count * 2) << " inner gap links) of "
            << total_links << " total links" << std::endl;
        }
        */

        // trim horizontal and vertical spaces to fit
        uint64_t min_x = width;
        uint64_t max_x = std::numeric_limits<uint64_t>::min(); // 0
        uint64_t min_y = height + path_space;
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

        const char *filename = args::get(png_out_file).c_str();
        png::encodeOneStep(filename, crop, crop_width, crop_height);

        return 0;
    }

    static Subcommand odgi_viz("viz", "Visualize a variation graph in 1D.",PIPELINE, 3, main_viz);

}
