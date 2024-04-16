#include "draw.hpp"
#include "split.hpp"

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

namespace algorithms {

void get_layout(const std::vector<double> &X,
                const std::vector<double> &Y,
                const PathHandleGraph &graph,
                const double& scale,
                const double& border,
                std::vector<std::vector<handle_t>>& weak_components,
                coord_range_2d_t& rendered_range,
                std::vector<coord_range_2d_t>& component_ranges) {

    // determine boundaries
    weak_components = algorithms::weakly_connected_component_vectors(&graph);
    double curr_y_offset = border;
    for (auto& component : weak_components) {
        component_ranges.emplace_back();
        auto& component_range = component_ranges.back();
        for (auto& handle : component) {
            uint64_t i = 2 * number_bool_packing::unpack_number(handle);
            for (uint64_t j = i; j <= i+1; ++j) {
                double x = X[j] * scale;
                double y = Y[j] * scale;
                component_range.include(x, y);
            }
        }
        component_range.x_offset = component_range.min_x - border;
        component_range.y_offset = curr_y_offset - component_range.min_y;
        curr_y_offset += component_range.height() + border;
    }


    // now examine the coordinates to determine our window size
    rendered_range = {0, 0, 0, 0};
    for (auto& component_range : component_ranges) {
        /*
        std::cerr << component_range.min_x << " "
                  << component_range.max_x << " "
                  << component_range.min_y << " "
                  << component_range.max_x << std::endl;
        */
        rendered_range.include(
            component_range.width() + 2 * border,
            component_range.max_y + component_range.y_offset + border);
    }

}

struct label_info_t {
    double x, y;
    std::string content;
    // Simple constructor for convenience
    label_info_t(double x, double y, std::string content) : x(x), y(y), content(std::move(content)) {}
};
bool is_too_close(double x, double y, const std::string& content, double threshold, std::vector<label_info_t>& placed_labels) {
    for (const auto& label : placed_labels) {
        if (label.content == content && std::abs(label.x - x) < threshold && std::abs(label.y - y) < threshold) {
            return true; // Found a label too close with the same content
        }
    }
    return false;
}

uint64_t node_hash(const nid_t& node_id) {
    uint64_t x = node_id;
    x = (~x) + (x << 21); // x = (x << 21) - x - 1;
    x = x ^ (x >> 24);
    x = (x + (x << 3)) + (x << 8); // x * 265
    x = x ^ (x >> 14);
    x = (x + (x << 2)) + (x << 4); // x * 21
    x = x ^ (x >> 28);
    x = x + (x << 31);
    return x;
}
bool keep_node(const nid_t& node_id, const float f) {
    // hash the node_id and check if it's accepted given our sparsification factor
    return node_hash(node_id) < std::numeric_limits<uint64_t>::max() * f;
}
// Define a struct to hold the coordinates for simplicity
struct Coordinates {
    double x1, y1, x2, y2;
};

// Function to adjust node length
Coordinates adjustNodeLength(double x1, double y1, double x2, double y2, double scale, double x_off, double y_off, double sparsification_factor) {
    // Apply scale and offsets to original coordinates
    x1 = (x1 * scale) - x_off;
    y1 = (y1 * scale) + y_off;
    x2 = (x2 * scale) - x_off;
    y2 = (y2 * scale) + y_off;

    // Calculate the original length
    double length = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));

    // Adjust length based on 1.0 / sparsification_factor
    double new_length = sparsification_factor == 0 ? length : length * (1.0 / sparsification_factor);

    // Calculate the midpoint
    double mid_x = (x1 + x2) / 2.0;
    double mid_y = (y1 + y2) / 2.0;

    // Calculate the unit vector for the direction
    double unit_x = (x2 - x1) / length;
    double unit_y = (y2 - y1) / length;

    // Calculate new endpoints using the new length
    double half_new_length = new_length / 2.0;
    double new_x1 = mid_x - half_new_length * unit_x;
    double new_y1 = mid_y - half_new_length * unit_y;
    double new_x2 = mid_x + half_new_length * unit_x;
    double new_y2 = mid_y + half_new_length * unit_y;

    // Return the new coordinates
    return Coordinates{new_x1, new_y1, new_x2, new_y2};
}
Coordinates adjustNodeEndpoints(const handle_t& handle, const std::vector<double>& X, const std::vector<double>& Y, double scale, double x_off, double y_off, double sparsification_factor, bool lengthen_left_nodes) {
    // Original coordinates
    uint64_t a = 2 * number_bool_packing::unpack_number(handle);
    double x1 = (X[a] * scale) - x_off;
    double y1 = (Y[a] * scale) + y_off;
    double x2 = (X[a + 1] * scale) - x_off;
    double y2 = (Y[a + 1] * scale) + y_off;

    // Calculate the original length
    double length = std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));

    // Adjust length based on 1.0 / sparsification_factor
    double new_length = !lengthen_left_nodes || sparsification_factor == 0 ? length : length * (1.0 / sparsification_factor);

    // Calculate the midpoint
    double mid_x = (x1 + x2) / 2.0;
    double mid_y = (y1 + y2) / 2.0;

    // Calculate the unit vector for the direction
    double unit_x = (x2 - x1) / length;
    double unit_y = (y2 - y1) / length;

    // Calculate new endpoints using the new length
    double half_new_length = new_length / 2.0;
    double new_x1 = mid_x - half_new_length * unit_x;
    double new_y1 = mid_y - half_new_length * unit_y;
    double new_x2 = mid_x + half_new_length * unit_x;
    double new_y2 = mid_y + half_new_length * unit_y;

    return Coordinates{new_x1, new_y1, new_x2, new_y2};
}

void draw_svg(std::ostream &out,
              const std::vector<double> &X,
              const std::vector<double> &Y,
              const PathHandleGraph &graph,
              const double& scale,
              const double& border,
			  const double& line_width,
			  std::vector<algorithms::color_t>& node_id_to_color,
              ska::flat_hash_map<handlegraph::nid_t, std::set<std::string>>& node_id_to_label_map,
              const float& sparsification_factor,
              const bool& lengthen_left_nodes) {

    std::vector<std::vector<handle_t>> weak_components;
    coord_range_2d_t rendered_range;
    std::vector<coord_range_2d_t> component_ranges;
    get_layout(X, Y, graph, scale, border, weak_components, rendered_range, component_ranges);

    double viewbox_x1 = rendered_range.min_x;
    double viewbox_x2 = rendered_range.max_x;
    double viewbox_y1 = rendered_range.min_y;
    double viewbox_y2 = rendered_range.max_y;
    
    double width = rendered_range.width();
    double height = rendered_range.height();

    std::vector<label_info_t> placed_labels;

    out << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    out << "<svg width=\"" << width << "\" height=\"" << height << "\" "
        << "viewBox=\"" << viewbox_x1 << " " << viewbox_y1
        << " " << width << " " << height << "\""
        << " xmlns=\"http://www.w3.org/2000/svg\">"
		// interferes with the styling of the lines
        //<< "<style type=\"text/css\">"
        //<< "line{stroke:black;stroke-width:1.0;stroke-opacity:1.0;stroke-linecap:round;};"
        //<< "</style>"
        << std::endl;

    auto range_itr = component_ranges.begin();
    for (auto& component : weak_components) {
        auto& range = *range_itr++;
        auto& x_off = range.x_offset;
        auto& y_off = range.y_offset;
		//const algorithms::color_t node_color = !node_id_to_color.empty() ? node_id_to_color[graph.get_id(handle)] : COLOR_BLACK;

        std::vector<handle_t> highlights;
        std::vector<handle_t> nodes_with_labels;

        for (auto& handle : component) {
			algorithms::color_t color = node_id_to_color.empty() ? COLOR_BLACK : node_id_to_color[graph.get_id(handle)];

            if (!(sparsification_factor == 0 || keep_node(graph.get_id(handle), sparsification_factor) || node_id_to_label_map.count(graph.get_id(handle)))) {
                continue; // Skip this node to output a lighter SVG (do not nodes with labels, if any)
            }

            Coordinates newEndpoints = adjustNodeEndpoints(handle, X, Y, scale, x_off, y_off, sparsification_factor, lengthen_left_nodes);

            if (color == COLOR_BLACK || color == COLOR_LIGHTGRAY) {
                out << "<line x1=\""
                    << newEndpoints.x1
                    << "\" x2=\""
                    << newEndpoints.x2
                    << "\" y1=\""
                    << newEndpoints.y1
                    << "\" y2=\""
                    << newEndpoints.y2
                    << "\" stroke=\"" << to_hexrgb(color)
                    << "\" stroke-width=\"" << line_width
                    << "\"/>"
                    << std::endl;
            } else {
                highlights.push_back(handle);
            }

            // Check if this is a node with a label
            if (node_id_to_label_map.count(graph.get_id(handle))){
                nodes_with_labels.push_back(handle);
            }
        }

        // Color highlights and put them after grey nodes to have colored nodes on top of grey ones
        for (auto& handle : highlights) {
            Coordinates newEndpoints = adjustNodeEndpoints(handle, X, Y, scale, x_off, y_off, sparsification_factor, lengthen_left_nodes);
            algorithms::color_t color = node_id_to_color.empty() ? COLOR_BLACK : node_id_to_color[graph.get_id(handle)];
            out << "<line x1=\""
                << newEndpoints.x1
                << "\" x2=\""
                << newEndpoints.x2
                << "\" y1=\""
                << newEndpoints.y1
                << "\" y2=\""
                << newEndpoints.y2
                << "\" stroke=\"" << to_hexrgb(color) //to_rgba(color) // with rgba, nodes are invisible with InkScape
                << "\" stroke-width=\"" << line_width
                << "\"/>"
                << std::endl;
        }

        // Render labels at the end, to have them on top of everything
        for (auto& handle : nodes_with_labels) {
            Coordinates newEndpoints = adjustNodeEndpoints(handle, X, Y, scale, x_off, y_off, sparsification_factor, lengthen_left_nodes);

            // Collect the labels that can be put without overlapping identical ones
            std::vector<std::string> labels;
            for (auto text : node_id_to_label_map[graph.get_id(handle)]){
                if (!is_too_close(newEndpoints.x2, newEndpoints.y2, text, 30.0, placed_labels)) {
                    labels.push_back(text);
                }
            }
            // Check if there is something to label
            if (!labels.empty()){
                out << "<text font-family=\"Arial\" font-size=\"20\" fill=\"#000000\" stroke=\"#000000\" y=\"" << newEndpoints.y2 << "\">";
                for (auto text : labels){
                    out << "<tspan x=\"" << newEndpoints.x2 << "\" dy=\"1.0em\">" << text << "</tspan>";
                    placed_labels.emplace_back(newEndpoints.x2, newEndpoints.y2, text); // Record the label's placement
                }
                out << "</text>"
                    << std::endl;
            }
        }
    }

    // todo, edges, paths, coverage, bins

    out << "</svg>" << std::endl;
}

std::vector<uint8_t> rasterize(const std::vector<double> &X,
                               const std::vector<double> &Y,
                               const PathHandleGraph &graph,
                               const double& scale,
                               const double& border,
                               uint64_t& width,
                               uint64_t& height,
                               const double& line_width,
                               const double& path_line_spacing,
                               bool color_paths,
                               std::vector<algorithms::color_t>& node_id_to_color) {

    std::vector<std::vector<handle_t>> weak_components;
    coord_range_2d_t rendered_range;
    std::vector<coord_range_2d_t> component_ranges;
    get_layout(X, Y, graph, scale, border, weak_components, rendered_range, component_ranges);

    double source_min_x = rendered_range.min_x;
    double source_min_y = rendered_range.min_y;
    double source_width = rendered_range.width();
    double source_height = rendered_range.height();

    std::vector<color_t> all_path_colors;
    if (color_paths) {
        graph.for_each_path_handle(
            [&](const path_handle_t& p) {
                all_path_colors.push_back(
                    hash_color(graph.get_path_name(p)));
                //std::cerr << graph.get_path_name(p) << " color " << all_path_colors.back() << std::endl;
            });
    }

    // determine height and width based on the width, if width = 0
    if (width == 0) {
        // Avoid too big images that would require too much memory
        width = std::min(100000.0, std::ceil(height * (source_width / source_height)));
    }
    //std::cerr << "source " << source_width << "×" << source_height << std::endl;
    //std::cerr << "raster " << width << "×" << height << std::endl;

    // an atomic image buffer that we can write into in parallel
    atomic_image_buf_t image(width, height,
                             source_width, source_height,
                             source_min_x, source_min_y);

    auto range_itr = component_ranges.begin();
	struct draw_target_t {
		xy_d_t xy0;
		xy_d_t xy1;
		algorithms::color_t color;
	};
	
    for (auto& component : weak_components) {
        auto& range = *range_itr++;
        auto& x_off = range.x_offset;
        auto& y_off = range.y_offset;
		std::vector<draw_target_t> highlights;
//#pragma omp parallel for
        for (uint64_t i = 0; i < component.size(); ++i) {
            const handle_t& handle = component[i];
            uint64_t a = 2 * number_bool_packing::unpack_number(handle);
            xy_d_t xy0 = {
                (X[a] * scale) - x_off,
                (Y[a] * scale) + y_off
            };
            xy0.into(source_min_x, source_min_y,
                     source_width, source_height,
                     2, 2,
                     width-4, height-4);
            xy_d_t xy1 = {
                (X[a + 1] * scale) - x_off,
                (Y[a + 1] * scale) + y_off
            };
            xy1.into(source_min_x, source_min_y,
                     source_width, source_height,
                     2, 2,
                     width-4, height-4);
            if (color_paths) {
                std::vector<color_t> path_colors;
                graph.for_each_step_on_handle(
                    handle,
                    [&](const step_handle_t& s) {
                        path_colors.push_back(
                            all_path_colors[as_integer(graph.get_path_handle_of_step(s))-1]);
                    });
                // for step on handle
                // get the path color
                wu_calc_rainbow(xy0, xy1, image, path_colors, path_line_spacing, line_width);
            } else {
                /*
                aaline(xy0, xy1,
                       COLOR_BLACK,
                       image,
                       line_width);
                */
                const algorithms::color_t node_color = !node_id_to_color.empty() ? node_id_to_color[graph.get_id(handle)] : COLOR_BLACK;

				// if gray or black color, otherwise save for later
				if (node_color == COLOR_BLACK || node_color == COLOR_LIGHTGRAY) {
					wu_calc_wide_line(xy0, xy1, node_color, image, line_width);
				} else {
					highlights.push_back({xy0, xy1, node_color});
				}
            }
        }
		// color highlights
		for (auto& highlight : highlights) {
			wu_calc_wide_line(highlight.xy0, highlight.xy1, highlight.color, image, line_width);
		}
    }

    // todo, edges, paths, coverage, bins
    
    return image.to_bytes();
}

void draw_png(const std::string& filename,
              const std::vector<double> &X,
              const std::vector<double> &Y,
              const PathHandleGraph &graph,
              const double& scale,
              const double& border,
              uint64_t width,
              uint64_t height,
              const double& line_width,
              const double& path_line_spacing,
              bool color_paths,
              std::vector<algorithms::color_t>& node_id_to_color) {
    auto bytes = rasterize(X, Y,
                           graph,
                           scale,
                           border,
                           width,
                           height,
                           line_width,
                           path_line_spacing,
                           color_paths,
                           node_id_to_color);
    png::encodeOneStep(filename.c_str(), bytes, width, height);
}

}
}
