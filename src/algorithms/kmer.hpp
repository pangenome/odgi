#ifndef DSGVG_KMER_HPP_INCLUDED
#define DSGVG_KMER_HPP_INCLUDED

#include <iostream>
#include <string>
#include <list>
#include "handle.hpp"
#include "position.hpp"

/** \file 
 * Functions for working with `kmers_t`'s in HandleGraphs.
 */

namespace dsgvg {

/// Stores a kmer in the context of a graph.
struct kmer_t {
    kmer_t(const std::string& s,
           const pos_t& b,
           const pos_t& e,
           const handle_t& c)
        : seq(s), begin(b), end(e), curr(c) { };
    /// the kmer
    std::string seq;
    /// our start position
    pos_t begin;
    /// Used in construction
    pos_t end; /// one past the (current) end of the kmer
    handle_t curr; /// the next handle we extend into
};

/// Print a kmer_t to a stream.
ostream& operator<<(ostream& out, const kmer_t& kmer);

namespace algorithms {

/// Iterate over all the kmers in the graph, running lambda on each
void for_each_kmer(const HandleGraph& graph, size_t k,
                   const function<void(const kmer_t&)>& lambda);

}

}

#endif
