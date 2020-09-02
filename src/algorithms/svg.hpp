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

namespace odgi {
namespace algorithms {

using namespace handlegraph;

void draw_svg(std::ostream &out,
              const std::vector<double> &X,
              const std::vector<double> &Y,
              const HandleGraph &graph,
              const double& scale,
              const double& border);

}
}
