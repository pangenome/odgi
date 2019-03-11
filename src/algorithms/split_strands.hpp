#ifndef VG_ALGORITHMS_SPLIT_STRANDS_HPP_INCLUDED
#define VG_ALGORITHMS_SPLIT_STRANDS_HPP_INCLUDED

/**
 * \file split_strands.hpp
 *
 * Defines algorithm for converting any graph into a single stranded graph.
 */

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/mutable_handle_graph.hpp>
#include "../utility.hpp"

#include <unordered_set>
#include <unordered_map>

namespace vg {
namespace algorithms {

using namespace std;
using namespace handlegraph;

    /// Fills a MutableHandleGraph 'into' with a graph that has the same sequence and path
    /// space as 'source', but all of the sequences are on the forward strand. This is
    /// accomplished by creating a new node for each node in the source graph with the reverse
    /// complement sequence. Returns a map that translates node IDs from 'into' to their
    /// node ID and orientation in 'source'. Reports an error and exits if 'into' is not
    /// empty.
    unordered_map<handlegraph::nid_t, pair<handlegraph::nid_t, bool>> split_strands(const HandleGraph* source,
                                                                                  MutableHandleGraph* into);

}
}

#endif
