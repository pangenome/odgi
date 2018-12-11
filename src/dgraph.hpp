//
//  dgraph.hpp
//  
//
//

#ifndef dgraph_hpp
#define dgraph_hpp

#include <cstdio>
#include <cstdint>
#include <vector>
#include <utility>
#include "sdsl/bit_vectors.hpp"
#include "sdsl/enc_vector.hpp"
#include "sdsl/dac_vector.hpp"
#include "sdsl/vlc_vector.hpp"
#include "sdsl/wavelet_trees.hpp"
#include "sdsl/csa_wt.hpp"
#include "sdsl/suffix_arrays.hpp"

namespace dankgraph {
    
    class SuccinctDynamicSequenceGraph {
        
    public:
        SuccinctDynamicSequenceGraph();
        ~SuccinctDynamicSequenceGraph();
        
        
    private:
        
        /// Encodes the topology of the graph.
        /// {ID, start edge list index, end edge list index, seq index}
        sdsl::int_vector<> graph_iv;
        
        /// Encodes a series of edges lists of nodes.
        /// {relative offset, orientation, next edge index}
        sdsl::int_vector<> edge_lists_iv;
        
        /// Encodes all of the sequences of all nodes and all paths in the graph.
        /// The node sequences occur in the same order as in graph_iv;
        sdsl::int_vector<> seq_iv;
        
        /// Same length as seq_iv. 1's indicate the beginning of a node's sequence.
        sdsl::bit_vector boundary_bv;
        
        /// Same length as seq_iv. 0's indicate that a base is still touched by some
        /// node or some path. 1's indicate that all nodes or paths that touch this
        /// base have been deleted.
        sdsl::bit_vector dead_bv;
        
        /// Encodes a self-balancing binary tree as integers. Consists of fixed-width
        /// records that have the following structure:
        /// {interval start, members index, parent index, left child index, right child index}
        /// Interval start variable indicates the start of a range in seq_iv (corresponding to
        /// a node, unless the node has been deleted), members index indicates the 1-based index
        /// of the first path membership record corresponding to this interval in
        /// path_membership_value_iv, and parent/left child/right child index indicates the
        /// topology of a binary tree search structure for these intervals. The indexes are 1-based
        /// with 0 indicating that the neighbor does not exist.
        sdsl::int_vector<> path_membership_range_iv;
        
        /// Encodes a series of linked lists. Consists of fixed-width records that have
        /// the following structure:
        /// {path id, rank, next index}
        /// Path ID indicates which path the node occurs on, rank indicates the ordinal
        /// position of this occurrence in the path, and next index indicates the 1-based
        /// index of the next occurrence of this node in this vector (or 0 if there is none)
        sdsl::int_vector<> path_membership_value_iv;
        
        /// Encodes the embedded paths of the graph. Each path is represented as a vector of
        /// sequence intervals. The sequence intervals are fixed-width records with the following
        /// structure:
        /// {inteval start, interval end}
        /// These values correspond to the 0-based indexes of an interval in seq_iv. The end index
        /// is past-the-last.
        /// The strand of this interval is given by the corresponding bit in the bit vector, with
        /// a 1 indicating reverse strand.
        std::vector<std::pair<sdsl::int_vector<>, sdsl::bit_vector>> paths;
        
        size_t dead_bases = 0;
        size_t deleted_nodes = 0;
        size_t deleted_edges = 0;
    };
}

#endif /* dgraph_hpp */
