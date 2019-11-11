#pragma once

/**
 * \file sparse_sort.hpp
 *
 * Defines a sparse matrix sort of the graph based on Mondriaan
 */

#include <iostream>
#include "../hash_map.hpp"
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>

#define GAINBUCKET_ARRAY
#define MONDRIAANVERSION "\"4.2.1\""
// hack to allow mondriaan to use variables named "new"
#define new mynew
extern "C" {
#include <Mondriaan.h>
}
#undef new

//#include "dynamic.hpp"
//#include "apply_bulk_modifications.hpp"
//#include "is_single_stranded.hpp"


namespace odgi {
namespace algorithms {

using namespace handlegraph;

void mondriaan_main(int argc, char** argv);

std::vector<handle_t> mondriaan_sort(const HandleGraph& graph);

}
}
