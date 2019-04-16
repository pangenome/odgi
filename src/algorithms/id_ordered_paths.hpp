#pragma once

#include <vector>
#include <utility>
#include <algorithm>
#include <limits>
#include <iostream>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/path_handle_graph.hpp>

namespace odgi {
namespace algorithms {

using namespace std;
using namespace handlegraph;

std::vector<path_handle_t> id_ordered_paths(const PathHandleGraph& g, bool rev = false);

}
}
