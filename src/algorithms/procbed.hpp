#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <random>
#include <fstream>
#include <omp.h>
#include "hash_map.hpp"
#include "progress.hpp"
#include "position.hpp"
#include "split.hpp"
#include "IITree.h"
#include <handlegraph/types.hpp>
#include <handlegraph/iteratee.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>

namespace odgi {

namespace algorithms {

using namespace handlegraph;

/// Subset and adjust the BED file to match the reference sub-ranges in the graph
void adjust_ranges(const PathHandleGraph& graph, const std::string& bed_targets);

}

}
