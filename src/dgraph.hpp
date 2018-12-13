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
#include <functional>
#include "path.hpp"
#include "handle_types.hpp"
#include "handle_helper.hpp"
#include "sdsl/bit_vectors.hpp"
//#include "sdsl/enc_vector.hpp"
//#include "sdsl/dac_vector.hpp"
//#include "sdsl/vlc_vector.hpp"
//#include "sdsl/wavelet_trees.hpp"
//#include "sdsl/csa_wt.hpp"
//#include "sdsl/suffix_arrays.hpp"

namespace dankgraph {


class SuccinctDynamicSequenceGraph {
        
public:
    SuccinctDynamicSequenceGraph();
    ~SuccinctDynamicSequenceGraph();
        
    /// Look up the handle for the node with the given ID in the given orientation
    handle_t get_handle(const id_t& node_id, bool is_reverse = false) const;
    
    /// Get the ID from a handle
    id_t get_id(const handle_t& handle) const;
    
    /// Get the orientation of a handle
    bool get_is_reverse(const handle_t& handle) const;
    
    /// Invert the orientation of a handle (potentially without getting its ID)
    handle_t flip(const handle_t& handle) const;
    
    /// Get the length of a node
    size_t get_length(const handle_t& handle) const;
    
    /// Get the sequence of a node, presented in the handle's local forward orientation.
    std::string get_sequence(const handle_t& handle) const;
    
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns true if we finished and false if we stopped early.
    bool follow_edges(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const;
    
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee
    /// returns false. Can be told to run in parallel, in which case stopping
    /// after a false return value is on a best-effort basis and iteration
    /// order is not defined.
    void for_each_handle(const std::function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
    
    /// Return the number of nodes in the graph
    /// TODO: can't be node_count because XG has a field named node_count.
    size_t node_size(void) const;
    
    /// Return the smallest ID in the graph, or some smaller number if the
    /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
    id_t min_node_id(void) const;
    
    /// Return the largest ID in the graph, or some larger number if the
    /// largest ID is unavailable. Return value is unspecified if the graph is empty.
    id_t max_node_id(void) const;
    
    
/**
 * This is the interface for a handle graph that stores embedded paths.
 */
    
//    ////////////////////////////////////////////////////////////////////////////
//    // Path handle interface that needs to be implemented
//    ////////////////////////////////////////////////////////////////////////////
//
//    /// Determine if a path name exists and is legal to get a path handle for.
//    bool has_path(const std::string& path_name) const;
//
//    /// Look up the path handle for the given path name.
//    /// The path with that name must exist.
//    path_handle_t get_path_handle(const std::string& path_name) const;
//
//    /// Look up the name of a path from a handle to it
//    std::string get_path_name(const path_handle_t& path_handle) const;
//
//    /// Returns the number of node occurrences in the path
//    size_t get_occurrence_count(const path_handle_t& path_handle) const;
//
//    /// Returns the number of paths stored in the graph
//    size_t get_path_count() const;
//
//    /// Execute a function on each path in the graph
//    // TODO: allow stopping early?
//    void for_each_path_handle(const std::function<void(const path_handle_t&)>& iteratee) const;
//
//    /// Get a node handle (node ID and orientation) from a handle to an occurrence on a path
//    handle_t get_occurrence(const occurrence_handle_t& occurrence_handle) const;
//
//    /// Get a handle to the first occurrence in a path.
//    /// The path MUST be nonempty.
//    occurrence_handle_t get_first_occurrence(const path_handle_t& path_handle) const;
//
//    /// Get a handle to the last occurrence in a path
//    /// The path MUST be nonempty.
//    occurrence_handle_t get_last_occurrence(const path_handle_t& path_handle) const;
//
//    /// Returns true if the occurrence is not the last occurence on the path, else false
//    bool has_next_occurrence(const occurrence_handle_t& occurrence_handle) const;
//
//    /// Returns true if the occurrence is not the first occurence on the path, else false
//    bool has_previous_occurrence(const occurrence_handle_t& occurrence_handle) const;
//
//    /// Returns a handle to the next occurrence on the path
//    occurrence_handle_t get_next_occurrence(const occurrence_handle_t& occurrence_handle) const;
//
//    /// Returns a handle to the previous occurrence on the path
//    occurrence_handle_t get_previous_occurrence(const occurrence_handle_t& occurrence_handle) const;
//
//    /// Returns a handle to the path that an occurrence is on
//    path_handle_t get_path_handle_of_occurrence(const occurrence_handle_t& occurrence_handle) const;
//
//    /// Returns the 0-based ordinal rank of a occurrence on a path
//    size_t get_ordinal_rank_of_occurrence(const occurrence_handle_t& occurrence_handle) const;
//
//    ////////////////////////////////////////////////////////////////////////////
//    // Additional optional interface with a default implementation
//    ////////////////////////////////////////////////////////////////////////////
//
//    /// Returns true if the given path is empty, and false otherwise
//    bool is_empty(const path_handle_t& path_handle) const;


/**
 * This is the interface for a handle graph that supports modification.
 */
    /*
     * Note: All operations may invalidate path handles and occurrence handles.
     */
    
    /// Create a new node with the given sequence and return the handle.
    handle_t create_handle(const std::string& sequence);

    /// Create a new node with the given id and sequence, then return the handle.
    handle_t create_handle(const std::string& sequence, const id_t& id);
    
    /// Remove the node belonging to the given handle and all of its edges.
    /// Does not update any stored paths.
    /// Invalidates the destroyed handle.
    /// May be called during serial for_each_handle iteration **ONLY** on the node being iterated.
    /// May **NOT** be called during parallel for_each_handle iteration.
    /// May **NOT** be called on the node from which edges are being followed during follow_edges.
    void destroy_handle(const handle_t& handle);
    
    /// Create an edge connecting the given handles in the given order and orientations.
    /// Ignores existing edges.
    void create_edge(const handle_t& left, const handle_t& right);
    
    /// Convenient wrapper for create_edge.
    inline void create_edge(const edge_t& edge) {
        create_edge(edge.first, edge.second);
    }
    
    /// Remove the edge connecting the given handles in the given order and orientations.
    /// Ignores nonexistent edges.
    /// Does not update any stored paths.
    void destroy_edge(const handle_t& left, const handle_t& right);
    
    /// Convenient wrapper for destroy_edge.
    inline void destroy_edge(const edge_t& edge) {
        destroy_edge(edge.first, edge.second);
    }
    
    /// Remove all nodes and edges. Does not update any stored paths.
    void clear(void);
    
    /// Swap the nodes corresponding to the given handles, in the ordering used
    /// by for_each_handle when looping over the graph. Other handles to the
    /// nodes being swapped must not be invalidated. If a swap is made while
    /// for_each_handle is running, it affects the order of the handles
    /// traversed during the current traversal (so swapping an already seen
    /// handle to a later handle's position will make the seen handle be visited
    /// again and the later handle not be visited at all).
    void swap_handles(const handle_t& a, const handle_t& b);
    
    /// Alter the node that the given handle corresponds to so the orientation
    /// indicated by the handle becomes the node's local forward orientation.
    /// Rewrites all edges pointing to the node and the node's sequence to
    /// reflect this. Invalidates all handles to the node (including the one
    /// passed). Returns a new, valid handle to the node in its new forward
    /// orientation. Note that it is possible for the node's ID to change.
    /// Does not update any stored paths. May change the ordering of the underlying
    /// graph.
    handle_t apply_orientation(const handle_t& handle);
    
    /// Split a handle's underlying node at the given offsets in the handle's
    /// orientation. Returns all of the handles to the parts. Other handles to
    /// the node being split may be invalidated. The split pieces stay in the
    /// same local forward orientation as the original node, but the returned
    /// handles come in the order and orientation appropriate for the handle
    /// passed in.
    /// Updates stored paths.
    std::vector<handle_t> divide_handle(const handle_t& handle, const std::vector<size_t>& offsets);

/**
 * This is the interface for a handle graph with embedded paths where the paths can be modified.
 * Note that if the *graph* can also be modified, the implementation will also
 * need to inherit from MutableHandleGraph, via the combination
 * MutablePathMutableHandleGraph interface.
 * TODO: This is a very limited interface at the moment. It will probably need to be extended.
 */
    
//    /**
//     * Destroy the given path. Invalidates handles to the path and its node occurrences.
//     */
//    void destroy_path(const path_handle_t& path);
//
//    /**
//     * Create a path with the given name. The caller must ensure that no path
//     * with the given name exists already, or the behavior is undefined.
//     * Returns a handle to the created empty path. Handles to other paths must
//     * remain valid.
//     */
//    path_handle_t create_path_handle(const std::string& name);
//
//    /**
//     * Append a visit to a node to the given path. Returns a handle to the new
//     * final occurrence on the path which is appended. Handles to prior
//     * occurrences on the path, and to other paths, must remain valid.
//     */
//    occurrence_handle_t append_occurrence(const path_handle_t& path, const handle_t& to_append);

/// These are the backing data structures that we use to fulfill the above functions

private:
    
    id_t max_id = 0;
    id_t min_id = std::numeric_limits<id_t>::max();
    
    inline uint8_t encode_nucleotide(const char& nt);
    inline char decode_nucleotide(const uint64_t& val);
    inline uint64_t complement_encoded_nucleotide(const uint64_t& val);
    inline size_t graph_iv_index(const handle_t& handle);
    inline uint64_t encode_edge_target(const handle_t& handle);
    inline handle decode_edge_target(const handle_t& handle);
    
    const static size_t GRAPH_RECORD_SIZE = 5;
    
    const static size_t GRAPH_ID_OFFSET = 0;
    const static size_t GRAPH_START_EDGES_OFFSET = 1;
    const static size_t GRAPH_END_EDGES_OFFSET = 2;
    const static size_t GRAPH_SEQ_START_OFFSET = 3;
    const static size_t GRAPH_SEQ_LENGTH_OFFSET = 4;
    
    const static size_t EDGE_RECORD_SIZE = 2;
    
    const static size_t EDGE_TRAV_OFFSET = 0;
    const static size_t EDGE_NEXT_OFFSET = 1;
    
    /// Encodes the topology of the graph.
    /// {ID, start edge list index, end edge list index, seq index, seq length}
    SuccinctDynamicVector graph_iv;

    /// Encodes a series of edges lists of nodes.
    /// {ID|orientation, next edge index}
    SuccinctDynamicVector edge_lists_iv;
    
    /// Encodes the 1-based offset of an ID in graph_iv in units of GRAPH_RECORD_SIZE.
    /// If no node with that ID exists, contains a 0.
    SuccinctDeque id_to_graph_iv;

    /// Encodes all of the sequences of all nodes and all paths in the graph.
    /// The node sequences occur in the same order as in graph_iv;
    SuccinctDynamicVector seq_iv;

    /// Same length as seq_iv. 1's indicate the beginning of a node's sequence.
    SuccinctDynamicVector boundary_bv;

    /// Same length as seq_iv. 0's indicate that a base is still touched by some
    /// node or some path. 1's indicate that all nodes or paths that touch this
    /// base have been deleted.
    SuccinctDynamicVector dead_bv;

    /// Encodes a self-balancing binary tree as integers. Consists of fixed-width
    /// records that have the following structure:
    /// {interval start, members index, parent index, left child index, right child index}
    /// Interval start variable indicates the start of a range in seq_iv (corresponding to
    /// a node, unless the node has been deleted), members index indicates the 1-based index
    /// of the first path membership record corresponding to this interval in
    /// path_membership_value_iv, and parent/left child/right child index indicates the
    /// topology of a binary tree search structure for these intervals. The indexes are 1-based
    /// with 0 indicating that the neighbor does not exist.
    SuccinctDynamicVector path_membership_range_iv;

    /// Encodes a series of linked lists. Consists of fixed-width records that have
    /// the following structure:
    /// {path id, rank, next index}
    /// Path ID indicates which path the node occurs on, rank indicates the ordinal
    /// position of this occurrence in the path, and next index indicates the 1-based
    /// index of the next occurrence of this node in this vector (or 0 if there is none)
    SuccinctDynamicVector path_membership_value_iv;

    /// Encodes the embedded paths of the graph. Each path is represented as three vectors
    /// starts, lengths, orientations
    /// The values in starts correspond to the 0-based indexes of an interval in seq_iv.
    /// The values in lengths are simply the length.
    /// The strand of this interval is given by the corresponding bit in orientations, with 1
    /// indicating reverse strand.
    std::vector<path_t> paths;

    size_t dead_bases;
    size_t deleted_nodes;
    size_t deleted_edges;
};
    
inline uint8_t SuccinctDynamicSequenceGraph::encode_nucleotide(const char& nt) {
    if (nt == 'a' || nt == 'A') {
        return 0;
    }
    else if (nt == 'c' || nt == 'C') {
        return 1;
    }
    else if (nt == 'g' || nt == 'G') {
        return 2;
    }
    else if (nt == 't' || nt == 'T') {
        return 3;
    }
    else {
        // all others, but probably N's
        return 4;
    }
}
    
inline uint64_t SuccinctDynamicSequenceGraph::complement_encoded_nucleotide(const uint64_t& val) {
    return val == 4 ? 4 : 3 - val;
}
    
inline uint8_t SuccinctDynamicSequenceGraph::decode_nucleotide(const uint64_t& val) {
    static const char* alphabet = "ACGTN";
    return alphabet[val];
}
    
inline size_t SuccinctDynamicSequenceGraph::graph_iv_index(const handle_t& handle) {
    return (id_to_graph_iv.get(get_id(handle)) - 1) * GRAPH_RECORD_SIZE;
}
    
inline uint64_t SuccinctDynamicSequenceGraph::encode_edge_target(const handle_t& handle) {
    return reinterpret_cast<uint64_t>(handle);
}
    
inline handle_t SuccinctDynamicSequenceGraph::decode_edge_target(const uint64_t& val) {
    return reinterpret_cast<handle_t>(val);
}
    
std::string reverse_complement(const std::string& seq) {
    std::string rev_comp(seq.size());
    for (size_t i = 0; i < seq.size(); i++) {
        char nt = seq.at(i);
        if (nt == 'a' || nt == 'A') {
            rev_comp[rev_comp.size() - i - 1] = 'T';
        }
        else if (nt == 'c' || nt == 'C') {
            rev_comp[rev_comp.size() - i - 1] = 'G';
        }
        else if (nt == 'g' || nt == 'G') {
            rev_comp[rev_comp.size() - i - 1] = 'C';
        }
        else if (nt == 't' || nt == 'T') {
            rev_comp[rev_comp.size() - i - 1] = 'A';
        }
        else {
            rev_comp[rev_comp.size() - i - 1] = 'N';
        }
    }
    return rev_comp;
}

} // end dankness

#endif /* dgraph_hpp */
