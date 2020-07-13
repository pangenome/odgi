#ifndef VG_ALGORITHMS_EXTEND_HPP_INCLUDED
#define VG_ALGORITHMS_EXTEND_HPP_INCLUDED

/**
 * \file extend.hpp
 *
 * Defines algorithm adding graph material from one graph into another
 */

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/mutable_handle_graph.hpp>
#include "utility.hpp"

namespace vg {
namespace algorithms {

using namespace std;
using namespace handlegraph;

    /// Adds the non-duplicative nodes and edges from 'source' to 'into'. Assumes that
    /// both graphs use the same ID space.
    void extend(const HandleGraph* source, MutableHandleGraph* into);

}
}

#endif
