#pragma once

#include <vector>
#include <utility>
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <iomanip> // std::setprecision
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include "progress.hpp"

namespace odgi {
    namespace algorithms {

        using namespace std;
        using namespace handlegraph;

        void bin_path_depth(const PathHandleGraph &graph,
                           const bool progress = false,
                           const uint64_t min_paths = 1,
                           const uint64_t min_depth = 1);
    }
}
