#include "svg.hpp"

namespace odgi {
namespace algorithms {

void draw_svg(std::ostream &out,
              const std::vector<double> &X,
              const std::vector<double> &Y,
              const HandleGraph &graph,
              const double& scale,
              const double& border) {
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::min();
    double max_y = std::numeric_limits<double>::min();

    // determine boundaries
    auto weak_components = algorithms::weakly_connected_component_vectors(&graph);
    std::vector<uint64_t> y_offsets;
    uint64_t y_offset = 0;
    for (auto& component : weak_components) {
        y_offsets.push_back(y_offset);
        for (auto& handle : component) {
            uint64_t i = number_bool_packing::unpack_number(handle);
            for (uint64_t j = i; j <= i+1; ++j) {
                double x = X[j] * scale;
                double y = Y[j] * scale + y_offset;
                if (x < min_x) min_x = x;
                if (x > max_x) max_x = x;
                if (y < min_y) min_y = y;
                if (y > max_y) max_y = y;
            }
        }
        //y_offset += max_y + border; // layout vertically
    }

    double viewbox_x1 = min_x;
    double viewbox_x2 = max_x;
    double viewbox_y1 = min_y;
    double viewbox_y2 = max_y;
    
    double width = viewbox_x2 - viewbox_x1;
    double height = viewbox_y2 - viewbox_y1;
    std::cerr << "width: " << width << std::endl;
    std::cerr << "height: " << height << std::endl;

    out << "<svg width=\"" << width << "\" height=\"" << height << "\" "
        << "viewBox=\"" << viewbox_x1 << " " << viewbox_y1
        << " " << width << " " << height << "\""
        << " xmlns=\"http://www.w3.org/2000/svg\">"
        << "<style type=\"text/css\">"
        << "line{stroke:black;stroke-width:1.0;stroke-opacity:1.0;stroke-linecap:round;}"
        << "</style>"
        << std::endl;

    auto y_offset_itr = y_offsets.begin();
    for (auto& component : weak_components) {
        uint64_t y_off = *y_offset_itr;
        ++y_offset_itr;
        for (auto& handle : component) {
            uint64_t a = 2 * number_bool_packing::unpack_number(handle);
            //std::cerr << a << ": " << X[a] << "," << Y[a] << " ------ " << X[a + 1] << "," << Y[a + 1] << std::endl;
            out << "<line x1=\""
                << X[a] * scale
                << "\" x2=\""
                << X[a + 1] * scale
                << "\" y1=\""
                << y_off + Y[a] * scale
                << "\" y2=\""
                << y_off + Y[a + 1] * scale
                << "\"/>"
                << std::endl;

        }
    }

    /*
    graph.for_each_handle(
        [&](const handle_t& handle) {
            uint64_t a = 2 * number_bool_packing::unpack_number(handle);
            //std::cerr << a << ": " << X[a] << "," << Y[a] << " ------ " << X[a + 1] << "," << Y[a + 1] << std::endl;
            out << "<line x1=\""
                << border + X[a] * scale
                << "\" x2=\""
                << border + X[a + 1] * scale
                << "\" y1=\""
                << border + Y[a] * scale
                << "\" y2=\""
                << border + Y[a + 1] * scale
                << "\"/>"
                << std::endl;
        });
    */

    // iterate through graph edges
    /*graph.for_each_edge([&](const edge_t &e) {
      uint64_t a = 2 * number_bool_packing::unpack_number(e.first);
      uint64_t b = 2 * number_bool_packing::unpack_number(e.second);

      //std::cerr << a << ": " << X[a] << "," << Y[a] << " ------ " << X[a + 1] << "," << Y[a + 1] << std::endl;
      out << "<line x1=\"" << X[a] * scale << "\" x2=\"" << X[a + 1] * scale
      << "\" y1=\"" << Y[a] * scale << "\" y2=\"" << Y[a + 1] * scale << "\"/>"
      << std::endl;

      //std::cerr << b << ": " << X[b] << "," << Y[b] << " ------ " << X[b + 1] << "," << Y[b + 1] << std::endl;
      out << "<line x1=\"" << X[b] * scale << "\" x2=\"" << X[b + 1] * scale
      << "\" y1=\"" << Y[b] * scale << "\" y2=\"" << Y[b + 1] * scale << "\"/>"
      << std::endl;
      });*/


    /* // to draw nodes
       for (uint64_t i = 0; i < n; ++i) {
       std::cout << "<circle cx=\"" << X[i*2]*scale << "\" cy=\"" << X[i*2+1]*scale << "\" r=\"1.0\"/>" << std::endl;
       }
    */

    out << "</svg>" << std::endl;
}

}
}
