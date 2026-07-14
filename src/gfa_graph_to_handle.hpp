#pragma once

/**
 * \file gfa_graph_to_handle.hpp
 *
 * Build a mutable handle graph from an already-decompressed GFAz graph.
 *
 * The heavy lifting of turning a .gfaz container into an in-memory GfaGraph is
 * done by GFAz (deserialize_compressed_data + decompress_gfa). This function
 * only translates that decoded GfaGraph into odgi's handle graph, mirroring
 * what gfa_to_handle() does for text GFA.
 */

#include "gfa_parser.hpp" // GfaGraph, NodeId and the columnar record types from GFAz
#include <cstdint>
#include <handlegraph/mutable_path_mutable_handle_graph.hpp>

namespace odgi {

/// Fill an (empty) handle graph from a decoded GFAz GfaGraph.
/// gfa_graph is consumed: its columns are released as they are used so the
/// decoded copy and the succinct graph do not both sit fully in memory.
void gfa_graph_to_handle(GfaGraph &gfa_graph,
                         handlegraph::MutablePathMutableHandleGraph *graph,
                         bool compact_ids,
                         uint64_t n_threads,
                         bool show_progress);

}
