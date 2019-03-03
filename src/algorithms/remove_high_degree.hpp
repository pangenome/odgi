#ifndef DSGVG_ALGORITHMS_REMOVE_HIGH_DEGREE_HPP_INCLUDED
#define DSGVG_ALGORITHMS_REMOVE_HIGH_DEGREE_HPP_INCLUDED

/**
 * \file remove_high_degree.hpp
 *
 * Defines a process that removes high-degree nodes from a graph
 */

#include <handlegraph/deletable_handle_graph.hpp>
#include <vector>
#include <iostream>

namespace odgi {
namespace algorithms {

using namespace handlegraph;

void remove_high_degree_nodes(DeletableHandleGraph& g, int max_degree);

}
}

#endif
