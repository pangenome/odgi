#pragma once

/**
 * \file sgd_layout.hpp
 *
 * Run SGD based graph layout
 */

#include <handlegraph/handle_graph.hpp>
#include <vector>
#include <random>
#include "weakly_connected_components.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

std::vector<double> sgd_layout(const HandleGraph& graph, uint64_t pivots, uint64_t t_max, double eps, double x_padding);

}
}



