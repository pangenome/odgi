#ifndef DSGVG_PRUNE_HPP_INCLUDED
#define DSGVG_PRUNE_HPP_INCLUDED

#include <iostream>
#include <vector>
#include <list>
#include <omp.h>
#include <handlegraph/util.hpp>
#include <handlegraph/handle_graph.hpp>
#include "position.hpp"
#include "odgi.hpp"

/** \file
 * Functions for working with `kmers_t`'s in HandleGraphs.
 */

namespace odgi {

using namespace handlegraph;

namespace algorithms {

/// How to treat paths damaged when prune removes nodes and/or edges.
enum class damaged_path_policy_t {
    split,          // break each path into subpaths at every removed node/edge
    drop_affected   // keep only fully-intact paths; drop the rest
};

/// One step of a recorded path in 12 bytes. The node id is kept full-width (nid_t) as two uint32
/// halves so the struct stays 4-byte aligned; only the length is packed (31 bits, guarded in
/// record_paths). Lengths let us name split subpaths in original path coordinates.
struct recorded_step_t {
    uint32_t id_low;
    uint32_t id_high;
    uint32_t length : 31;
    uint32_t is_rev : 1;
    nid_t node_id() const { return (nid_t)(((uint64_t)id_high << 32) | id_low); }
    void set_node_id(nid_t v) { id_low = (uint32_t)v; id_high = (uint32_t)((uint64_t)v >> 32); }
};

/// A path snapshot taken while the graph is still intact.
struct recorded_path_t {
    std::string name;
    bool is_circular = false;
    std::vector<recorded_step_t> steps;
};

struct prune_path_rebuild_stats_t {
    uint64_t paths_intact = 0;     // re-embedded whole (untouched by pruning)
    uint64_t paths_split = 0;      // affected paths that yielded >=1 subpath
    uint64_t subpaths_created = 0; // total subpaths from split paths
    uint64_t paths_dropped = 0;    // nothing left, or dropped by drop_affected
};

/// Snapshot every path so it can be rebuilt after pruning. Must run before any destroy; we
/// snapshot rather than graph_t::copy because that copy corrupts path linkage.
std::vector<recorded_path_t> record_paths(const graph_t& graph);

/// Re-embed recorded paths into the already-pruned `subgraph` (which must have no paths). A run
/// breaks at any removed node or removed edge. split emits every maximal run (whole-path runs keep
/// the name, broken pieces are PATH:START-END); drop_affected keeps a path only if fully intact.
prune_path_rebuild_stats_t rebuild_pruned_paths(const std::vector<recorded_path_t>& recorded_paths,
                                                graph_t& subgraph,
                                                const damaged_path_policy_t policy);

/// Record a <=k-length walk in the context of a graph.
struct walk_t {
    walk_t(uint16_t l,
           const pos_t& b,
           const pos_t& e,
           const handle_t& c,
           uint16_t f)
        : length(l), begin(b), end(e), curr(c), forks(f) { };
    /// our start position
    pos_t begin;
    pos_t end; /// one past the (current) end of the kmer
    handle_t curr; /// the next handle we extend into
    uint16_t forks; /// how many branching edge crossings we took to get here
    uint16_t length; /// how far we've been
};

/// Iterate over all the walks up to length k, adding edges which 
std::vector<edge_t> find_edges_to_prune(const HandleGraph& graph, size_t k, size_t edge_max, int n_threads);

}

}

#endif
