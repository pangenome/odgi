#include "draw.hpp"

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
        component_range.y_offset = curr_y_offset -component_range.min_y;
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


void draw_svg(std::ostream &out,
              const std::vector<double> &X,
              const std::vector<double> &Y,
              const PathHandleGraph &graph,
              const double& scale,
              const double& border) {

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

    out << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    out << "<svg width=\"" << width << "\" height=\"" << height << "\" "
        << "viewBox=\"" << viewbox_x1 << " " << viewbox_y1
        << " " << width << " " << height << "\""
        << " xmlns=\"http://www.w3.org/2000/svg\">"
        << "<style type=\"text/css\">"
        << "line{stroke:black;stroke-width:1.0;stroke-opacity:1.0;stroke-linecap:round;};"
        << "</style>"
        << std::endl;

    auto range_itr = component_ranges.begin();
    for (auto& component : weak_components) {
        auto& range = *range_itr++;
        auto& x_off = range.x_offset;
        auto& y_off = range.y_offset;
        for (auto& handle : component) {
            uint64_t a = 2 * number_bool_packing::unpack_number(handle);
            out << "<line x1=\""
                << (X[a] * scale) - x_off
                << "\" x2=\""
                << (X[a + 1] * scale) - x_off
                << "\" y1=\""
                << (Y[a] * scale) + y_off
                << "\" y2=\""
                << (Y[a + 1] * scale) + y_off
                << "\"/>"
                << std::endl;

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
                               const double& line_overlay,
                               const double& path_line_spacing,
                               bool color_paths) {

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
            });
    }

    // determine height and width based on the width, if width = 0
    if (width == 0) {
        width = std::ceil(height * (source_width / source_height));
    }
    //std::cerr << "source " << source_width << "×" << source_height << std::endl;
    //std::cerr << "raster " << width << "×" << height << std::endl;

    // an atomic image buffer that we can write into in parallel
    atomic_image_buf_t image(width, height,
                             source_width, source_height,
                             source_min_x, source_min_y);

    auto range_itr = component_ranges.begin();
    for (auto& component : weak_components) {
        auto& range = *range_itr++;
        auto& x_off = range.x_offset;
        auto& y_off = range.y_offset;
#pragma omp parallel for
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
                            all_path_colors[as_integer(graph.get_path_handle_of_step(s))]);
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
                wu_calc_wide_line(xy0, xy1, COLOR_BLACK, image, line_width);
            }
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
              const double& line_overlay,
              const double& path_line_spacing,
              bool color_paths) {
    auto bytes = rasterize(X, Y,
                           graph,
                           scale,
                           border,
                           width,
                           height,
                           line_width,
                           line_overlay,
                           path_line_spacing,
                           color_paths);
    png::encodeOneStep(filename.c_str(), bytes, width, height);
}

}
}
