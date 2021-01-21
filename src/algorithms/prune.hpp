#ifndef DSGVG_PRUNE_HPP_INCLUDED
#define DSGVG_PRUNE_HPP_INCLUDED

#include <iostream>
#include <vector>
#include <list>
#include <omp.h>
#include <handlegraph/util.hpp>
#include <handlegraph/handle_graph.hpp>
#include "position.hpp"

/** \file 
 * Functions for working with `kmers_t`'s in HandleGraphs.
 */

namespace odgi {

using namespace handlegraph;

namespace algorithms {

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
