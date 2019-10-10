//
//  odgi
//
//  graph.hpp
//
//  main dynamic compact graph definition
//

#pragma once

#include <cstdio>
#include <cstdint>
#include <vector>
#include <utility>
#include <functional>
#include <handlegraph/types.hpp>
#include <handlegraph/iteratee.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/mutable_handle_graph.hpp>
#include <handlegraph/mutable_path_handle_graph.hpp>
#include <handlegraph/mutable_path_mutable_handle_graph.hpp>
#include <handlegraph/deletable_handle_graph.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
#include "dynamic.hpp"
#include "dynamic_types.hpp"
#include "dna.hpp"
#include "hash_map.hpp"
#include "node.hpp"

namespace odgi {

using namespace handlegraph;

// Resolve ambiguous nid_t typedef by putting it in our namespace.
using nid_t = handlegraph::nid_t;

class graph_t : public MutablePathDeletableHandleGraph {

public:

    graph_t(void) {
        // set up initial delimiters
        deleted_node_bv.push_back(1);
    }

    ~graph_t(void) { clear(); }

    /// Method to check if a node exists by ID
    bool has_node(nid_t node_id) const;
    
    /// Look up the handle for the node with the given ID in the given orientation
    handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const;

    /// Get the ID from a handle
    nid_t get_id(const handle_t& handle) const;
    
    /// Get the orientation of a handle
    bool get_is_reverse(const handle_t& handle) const;
    
    /// Invert the orientation of a handle (potentially without getting its ID)
    handle_t flip(const handle_t& handle) const;
    
    /// Get the length of a node
    size_t get_length(const handle_t& handle) const;
    
    /// Get the sequence of a node, presented in the handle's local forward orientation.
    std::string get_sequence(const handle_t& handle) const;
    
protected:
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns true if we finished and false if we stopped early.
    bool follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const;
    
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee
    /// returns false. Can be told to run in parallel, in which case stopping
    /// after a false return value is on a best-effort basis and iteration
    /// order is not defined.
    bool for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
    
public:

    /// Return the number of nodes in the graph
    /// TODO: can't be node_count because XG has a field named node_count.
    size_t get_node_count(void) const;
    
    /// Return the smallest ID in the graph, or some smaller number if the
    /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
    nid_t min_node_id(void) const;
    
    /// Return the largest ID in the graph, or some larger number if the
    /// largest ID is unavailable. Return value is unspecified if the graph is empty.
    nid_t max_node_id(void) const;

    /// Set a minimum id to increment the id space by, used as a hint during construction.
    /// May have no effect on a backing implementation.
    void set_id_increment(const nid_t& min_id);
    
    /// Get a handle from a Visit Protobuf object.
    /// Must be using'd to avoid shadowing.
    //handle_t get_handle(const Visit& visit) const;
    
    ////////////////////////////////////////////////////////////////////////////
    // Additional optional interface with a default implementation
    ////////////////////////////////////////////////////////////////////////////
    
    /// Get the number of edges on the right (go_left = false) or left (go_left
    /// = true) side of the given handle. The default implementation is O(n) in
    /// the number of edges returned, but graph implementations that track this
    /// information more efficiently can override this method.
    size_t get_degree(const handle_t& handle, bool go_left) const;
    
    /// Get the locally forward version of a handle
    handle_t forward(const handle_t& handle) const;
    
    /// A pair of handles can be used as an edge. When so used, the handles have a
    /// canonical order and orientation.
    edge_t edge_handle(const handle_t& left, const handle_t& right) const;
    
    /// Such a pair can be viewed from either inward end handle and produce the
    /// outward handle you would arrive at.
    handle_t traverse_edge_handle(const edge_t& edge, const handle_t& left) const;
    
    ////////////////////////////////////////////////////////////////////////////
    // Path handle interface
    ////////////////////////////////////////////////////////////////////////////
    
    /// Determine if a path name exists and is legal to get a path handle for.
    bool has_path(const std::string& path_name) const;
    
    /// Look up the path handle for the given path name.
    /// The path with that name must exist.
    path_handle_t get_path_handle(const std::string& path_name) const;
    
    /// Look up the name of a path from a handle to it
    std::string get_path_name(const path_handle_t& path_handle) const;
    
    /// Returns the number of node steps in the path
    size_t get_step_count(const path_handle_t& path_handle) const;

    /// Returns the number of paths stored in the graph
    size_t get_path_count(void) const;

    /// Returns a vector of all steps of a node on paths. Optionally restricts to
    /// steps that match the handle in orientation.
    std::vector<step_handle_t> steps_of_handle(const handle_t& handle,
                                               bool match_orientation = false) const;
    
protected:
    
    /// Execute a function on each path in the graph
    bool for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const;

    /// Enumerate the path steps on a given handle (strand agnostic)
    bool for_each_step_on_handle_impl(const handle_t& handle, const std::function<bool(const step_handle_t&)>& iteratee) const;
    
public:

    /// Returns the number of node steps on the handle
    size_t get_step_count(const handle_t& handle) const;
    
    /// Get a node handle (node ID and orientation) from a handle to an step on a path
    handle_t get_handle_of_step(const step_handle_t& step_handle) const;

    /// Get a path handle (path ID) from a handle to an step on a path
    path_handle_t get_path(const step_handle_t& step_handle) const;
    
    /// Get a handle to the first step in a path.
    /// The path MUST be nonempty.
    step_handle_t path_begin(const path_handle_t& path_handle) const;

    /// Get a handle to a fictitious handle one past the end of the path
    step_handle_t path_end(const path_handle_t& path_handle) const;

    /// Get a handle to the last step, which is arbitrary in the case of a circular path
    step_handle_t path_back(const path_handle_t& path_handle) const;

    /// Get a handle to a fictitious handle one past the start of the path
    step_handle_t path_front_end(const path_handle_t& path_handle) const;

    /// Returns true if the step handle is a front end magic handle
    bool is_path_front_end(const step_handle_t& step_handle) const;

    /// Returns true if the step handle is an end magic handle
    bool is_path_end(const step_handle_t& step_handle) const;
    
    /// Returns true if the step is not the last step on the path, else false
    bool has_next_step(const step_handle_t& step_handle) const;
    
    /// Returns true if the step is not the first step on the path, else false
    bool has_previous_step(const step_handle_t& step_handle) const;
    
    /// Returns a handle to the next step on the path
    step_handle_t get_next_step(const step_handle_t& step_handle) const;
    
    /// Returns a handle to the previous step on the path
    step_handle_t get_previous_step(const step_handle_t& step_handle) const;
    
    /// Returns a handle to the path that an step is on
    path_handle_t get_path_handle_of_step(const step_handle_t& step_handle) const;
    
    /// Returns the 0-based ordinal rank of a step on a path
    size_t get_ordinal_rank_of_step(const step_handle_t& step_handle) const;

    /// Returns true if the given path is empty, and false otherwise
    bool is_empty(const path_handle_t& path_handle) const;

    /// Loop over all the steps along a path, from first through last
    void for_each_step_in_path(const path_handle_t& path, const std::function<void(const step_handle_t&)>& iteratee) const;

    /// Returns true if the path is circular
    bool get_is_circular(const path_handle_t& path_handle) const;

    /// Set if the path is circular or not
    void set_circularity(const path_handle_t& path_handle, bool circular);
    
    /// Create a new node with the given sequence and return the handle.
    handle_t create_handle(const std::string& sequence);

    /// Create a new node with the given id and sequence, then return the handle.
    handle_t create_handle(const std::string& sequence, const nid_t& id);

    /// Create a "hidden" node which might carry parts of paths that traversed deleted portions of the graph
    handle_t create_hidden_handle(const std::string& sequence);

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

    /// Check if an edge exists
    bool has_edge(const handle_t& left, const handle_t& right) const;

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

    /// Remove all stored paths
    void clear_paths(void);
    
    /// Swap the nodes corresponding to the given handles, in the ordering used
    /// by for_each_handle when looping over the graph. Other handles to the
    /// nodes being swapped must not be invalidated. If a swap is made while
    /// for_each_handle is running, it affects the order of the handles
    /// traversed during the current traversal (so swapping an already seen
    /// handle to a later handle's position will make the seen handle be visited
    /// again and the later handle not be visited at all).
    void swap_handles(const handle_t& a, const handle_t& b);

    /// Reorder the graph's internal structure to match that given.
    /// Optionally compact the id space of the graph to match the ordering, from 1->|ordering|.
    void apply_ordering(const std::vector<handle_t>& order, bool compact_ids = false);

    /// Organize the graph for better performance and memory use
    void optimize(bool allow_id_reassignment = true);

    /// Reorder the graph's paths as given.
    void apply_path_ordering(const std::vector<path_handle_t>& order);
    
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
    
    /// Specialization of divide_handle for a single division point
    inline std::pair<handle_t, handle_t> divide_handle(const handle_t& handle, size_t offset) {
        auto parts = divide_handle(handle, std::vector<size_t>{offset});
        return std::make_pair(parts.front(), parts.back());
    }

    handle_t combine_handles(const std::vector<handle_t>& handles);

/**
 * This is the interface for a handle graph with embedded paths where the paths can be modified.
 * Note that if the *graph* can also be modified, the implementation will also
 * need to inherit from MutableHandleGraph, via the combination
 * MutablePathMutableHandleGraph interface.
 * TODO: This is a very limited interface at the moment. It will probably need to be extended.
 */
    
    /**
     * Destroy the given path. Invalidates handles to the path and its node steps.
     */
    void destroy_path(const path_handle_t& path);

    /**
     * Create a path with the given name. The caller must ensure that no path
     * with the given name exists already, or the behavior is undefined.
     * Returns a handle to the created empty path. Handles to other paths must
     * remain valid.
     */
    path_handle_t create_path_handle(const std::string& name,
                                     bool is_circular = false);

    /**
     * Append a visit to a node to the given path. Returns a handle to the new
     * final step on the path which is appended. Handles to prior
     * steps on the path, and to other paths, must remain valid.
     */
    step_handle_t prepend_step(const path_handle_t& path, const handle_t& to_append);
    
    /**
     * Append a visit to a node to the given path. Returns a handle to the new
     * final step on the path which is appended. Handles to prior
     * steps on the path, and to other paths, must remain valid.
     */
    step_handle_t append_step(const path_handle_t& path, const handle_t& to_append);

    /**
     * Insert a visit to a node to the given path between the given steps.
     * Returns a handle to the new step on the path which is appended.
     * Handles to prior steps on the path, and to other paths, must remain valid.
     */
    step_handle_t insert_step(const step_handle_t& before, const step_handle_t& after, const handle_t& to_insert);
    
    /// Set the step to the given handle, possibly re-linking and cleaning up if needed
    step_handle_t set_step(const step_handle_t& step_handle, const handle_t& handle);

    /// Replace the path range with the new segment
    std::pair<step_handle_t, step_handle_t> rewrite_segment(const step_handle_t& segment_begin,
                                                            const step_handle_t& segment_end,
                                                            const std::vector<handle_t>& new_segment);

    /// A helper function to visualize the state of the graph
    void display(void) const;

    /// Convert to GFA (for debugging)
    void to_gfa(std::ostream& out) const;

    /// Serialize
    uint64_t serialize(std::ostream& out);

    /// Load
    void load(std::istream& in);

/// These are the backing data structures that we use to fulfill the above functions

private:

    /// Records the handle to node_id mapping
    /// Use the special value "0" to indicate deleted nodes so that
    /// handle references in the id_map and elsewhere are not immediately destroyed
    std::vector<node_t> node_v;
    /// Mark deleted nodes here for translating graph ids into internal ranks
    suc_bv deleted_node_bv;
    uint64_t _deleted_node_count = 0;
    /// efficient id to handle/sequence conversion
    nid_t _max_node_id = 0;
    nid_t _min_node_id = 0;
    nid_t _id_increment = 0;
    /// records nodes that are hidden, but used to compactly store path sequence that has been removed from the node space
    hash_set<uint64_t> graph_id_hidden_set;

    /// edge type conversion
    /// 1 = fwd->fwd, 2 = fwd->rev, 3 = rev->fwd, 4 = rev->rev
    struct edge_helper {
        inline static uint8_t pack(bool on_rev, bool other_rev, bool to_curr) {
            return on_rev | (other_rev << 1) | (to_curr << 2);
        }
        inline static uint8_t unpack_on_rev(uint8_t edge) {
            return edge & 1;
        }
        inline static uint8_t unpack_other_rev(uint8_t edge) {
            return edge & (1 << 1);
        }
        inline static uint8_t unpack_to_curr(uint8_t edge) {
            return edge & (1 << 2);
        }
    };

    inline void canonicalize_edge(handle_t& left, handle_t& right) const {
        if (number_bool_packing::unpack_bit(left) && number_bool_packing::unpack_bit(right)
            || ((number_bool_packing::unpack_bit(left) || number_bool_packing::unpack_bit(right)) && as_integer(left) > as_integer(right))) {
            std::swap(left, right);
            left = number_bool_packing::toggle_bit(left);
            right = number_bool_packing::toggle_bit(right);
        }
    }
    
    struct path_metadata_t {
        uint64_t length;
        step_handle_t first;
        step_handle_t last;
        std::string name;
        bool is_circular = false;
    };
    /// maps between path identifier and the start, end, and length of the path
    std::vector<path_metadata_t> path_metadata_v;

    /// Links path names to handles
    string_hash_map<std::string, uint64_t> path_name_map;

    /// A helper to record the number of live nodes
    uint64_t _node_count = 0;

    /// A counter that records the number of hidden nodes
    uint64_t _hidden_count = 0;

    /// A helper to record the number of live edges
    uint64_t _edge_count = 0;

    /// A helper to record the number of live paths
    uint64_t _path_count = 0;

    /// A helper to record the next path handle (path deletions are hard because of our path FM-index)
    uint64_t _path_handle_next = 0;

    /// Helper to convert between edge storage and actual id
    uint64_t edge_delta_to_id(uint64_t left, uint64_t delta) const;

    /// Helper to convert between ids and stored edge
    uint64_t edge_to_delta(const handle_t& left, const handle_t& right) const;

    /// Helper to simplify removal of path handle records
    void destroy_path_handle_records(uint64_t i);

    /// Helper to create the internal records for the step
    step_handle_t create_step(const path_handle_t& path, const handle_t& handle);

    /// Helper to destroy the internal records for the step
    void destroy_step(const step_handle_t& step_handle);

    /// Helper to stitch up partially built paths
    void link_steps(const step_handle_t& from, const step_handle_t& to);

    /// Decrement the step rank references for this step
    void decrement_rank(const step_handle_t& step_handle);

    /// Compact away the deleted nodes info
    //void rebuild_id_handle_mapping(void);

    /// Set the handle sequence
    void set_handle_sequence(const handle_t& handle, const std::string& seq);

    /// get the backing node rank for a given node id
    uint64_t get_node_rank(const nid_t& node_id) const;

};

const static uint64_t path_begin_marker = 0;
const static uint64_t path_end_marker = 1;

} // end dankness
