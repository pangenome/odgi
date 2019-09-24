#ifndef DG_GFA_TO_HANDLE_HPP_INCLUDED
#define DG_GFA_TO_HANDLE_HPP_INCLUDED

/**
 * \file gfa_to_handle.hpp
 *
 * Contains a method to construct a mutable handle graph out of a GFA file
 *
 */

#include "gfakluge.hpp"
#include <iostream>
#include <limits>
#include <handlegraph/mutable_path_mutable_handle_graph.hpp>

namespace odgi {

    /// Fills a handle graph with an instantiation of a sequence graph from a GFA file.
    /// Handle graph must be empty when passed into function.
    void gfa_to_handle(const string& gfa_filename,
                       handlegraph::MutablePathMutableHandleGraph* graph,
                       bool show_progress = false);

}
#endif
