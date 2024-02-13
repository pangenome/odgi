#pragma once

#include <vector>
#include <string>
#include <cassert>
#include <algorithm>
#include <random>
#include <set>
#include <thread>
#include <atomic>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/handle_graph.hpp>
#include "weakly_connected_components.hpp"
#include <sdsl/bit_vectors.hpp>
#include <handlegraph/util.hpp>
#include "lodepng.h"
#include "atomic_image.hpp"
//#include "layout.hpp" // for callback interaction with succinct Layout

namespace odgi {

namespace png {

void encodeOneStep(const char *filename, std::vector<unsigned char> &image, unsigned width, unsigned height);
void encodeTwoSteps(const char *filename, std::vector<unsigned char> &image, unsigned width, unsigned height);
void encodeWithState(const char *filename, std::vector<unsigned char> &image, unsigned width, unsigned height);

}

namespace algorithms {

using namespace handlegraph;


struct coord_range_2d_t {
    double min_x = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::min();
    double min_y = std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::min();
    double x_offset = 0;
    double y_offset = 0;
    double width() {
        return max_x - min_x;
    }
    double height() {
        return max_y - min_y;
    }
    void include(const double& x,
                 const double& y) {
        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
    }
};


void get_layout(const std::vector<double> &X,
                const std::vector<double> &Y,
                const PathHandleGraph &graph,
                const double& scale,
                const double& border,
                std::vector<std::vector<handle_t>>& weak_components,
                coord_range_2d_t& rendered_range,
                std::vector<coord_range_2d_t>& component_ranges);

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
              const bool& lengthen_left_nodes);

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
                               std::vector<algorithms::color_t>& node_id_to_color);

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
              std::vector<algorithms::color_t>& node_id_to_color);



}

}
