#include "svg.hpp"

namespace odgi {
namespace algorithms {

void draw_svg(std::ostream &out,
              const std::vector<double> &X,
              const std::vector<double> &Y,
              const HandleGraph &graph,
              const double& scale,
              const double& border) {

    // determine boundaries
    auto weak_components = algorithms::weakly_connected_component_vectors(&graph);
    std::vector<coord_range_2d_t> component_ranges;
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
    coord_range_2d_t rendered_range = {0, 0, 0, 0};
    for (auto& component_range : component_ranges) {
        std::cerr << component_range.min_x << " "
                  << component_range.max_x << " "
                  << component_range.min_y << " "
                  << component_range.max_x << std::endl;
        rendered_range.include(
            component_range.width() + 2 * border,
            component_range.max_y + component_range.y_offset + border);
    }

    

    double viewbox_x1 = rendered_range.min_x;
    double viewbox_x2 = rendered_range.max_x;
    double viewbox_y1 = rendered_range.min_y;
    double viewbox_y2 = rendered_range.max_y;
    
    double width = rendered_range.width();
    double height = rendered_range.height();

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
        uint64_t x_off = range.x_offset;
        uint64_t y_off = range.y_offset;
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

}
}
