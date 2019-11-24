#pragma once

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/mutable_handle_graph.hpp>
#include "hash_map.hpp"
#include "topological_sort.hpp"
#include "split_strands.hpp"
#include "dagify.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

std::vector<handle_t> dagify_sort(const HandleGraph& base, MutableHandleGraph& split, MutableHandleGraph& into);

}
}
