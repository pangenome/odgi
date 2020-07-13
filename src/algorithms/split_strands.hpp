#pragma once

/**
 * \file split_strands.hpp
 *
 * Defines algorithm for converting any graph into a single stranded graph.
 */

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/mutable_handle_graph.hpp>
//#include "utility.hpp"
#include "dna.hpp"
#include "hash_map.hpp"
#include <utility>
#include <unordered_set>
#include <unordered_map>

namespace odgi {
namespace algorithms {

using namespace handlegraph;

    /// Fills a MutableHandleGraph 'into' with a graph that has the same sequence and path
    /// space as 'source', but all of the sequences are on the forward strand. This is
    /// accomplished by creating a new node for each node in the source graph with the reverse
    /// complement sequence. Returns a map that translates node IDs from 'into' to their
    /// node ID and orientation in 'source'. Reports an error and exits if 'into' is not
    /// empty.
ska::flat_hash_map<handlegraph::nid_t, std::pair<handlegraph::nid_t, bool>> split_strands(const HandleGraph* source,
                                                                                          MutableHandleGraph* into);

}
}
