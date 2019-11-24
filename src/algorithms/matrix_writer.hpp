#pragma once

/**
 * \file matrix_writer.hpp
 *
 * Defines a function to write a graph as a sparse matrix
 */

#include <iostream>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/util.hpp>

namespace odgi {
namespace algorithms {

using namespace handlegraph;

void write_as_sparse_matrix(std::ostream& out, const PathHandleGraph& graph, bool weight_by_edge_depth, bool weight_by_edge_delta);

}
}
