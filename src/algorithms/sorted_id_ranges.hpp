#ifndef VG_ALGORITHMS_SORTED_ID_RANGES_HPP_INCLUDED
#define VG_ALGORITHMS_SORTED_ID_RANGES_HPP_INCLUDED

/**
 * \file sorted_id_ranges.hpp
 *
 * Defines an algorithm to get the ranges of IDs covered by a HandleGraph.
 */

#include <vector>
#include <utility>

#include <handlegraph/handle_graph.hpp>

namespace vg {
namespace algorithms {

using namespace std;
using namespace handlegraph;

/// Get a sorted list of inclusive ranges of IDs used in the given HandleGraph.
vector<pair<handlegraph::nid_t, handlegraph::nid_t>> sorted_id_ranges(const HandleGraph* graph);


}
}

#endif
