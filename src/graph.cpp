//
//  graph.cpp
//  

#include "graph.hpp"

namespace dankgraph {

graph_t::graph_t(void) { }
graph_t::~graph_t(void) { }

/// Look up the handle for the node with the given ID in the given orientation
handle_t graph_t::get_handle(const id_t& node_id, bool is_reverse) const {
    return handle_helper::pack(graph_id_wt.select(1, node_id), is_reverse);
}
    
/// Get the ID from a handle
id_t graph_t::get_id(const handle_t& handle) const {
    return graph_id_wt.at(handle_helper::unpack_number(handle));
}
    
/// Get the orientation of a handle
bool graph_t::get_is_reverse(const handle_t& handle) const {
    return handle_helper::unpack_bit(handle);
}
    
/// Invert the orientation of a handle (potentially without getting its ID)
handle_t graph_t::flip(const handle_t& handle) const {
    return handle_helper::toggle_bit(handle);
}
    
/// Get the length of a node
size_t graph_t::get_length(const handle_t& handle) const {
    uint64_t offset = handle_helper::unpack_number(handle);
    return boundary_bv.select(offset+1) - boundary_bv.select(offset);
}
    
/// Get the sequence of a node, presented in the handle's local forward orientation.
std::string graph_t::get_sequence(const handle_t& handle) const {
    std::string seq;
    uint64_t offset = handle_helper::unpack_number(handle);
    uint64_t start = boundary_bv.select(offset);
    uint64_t end = boundary_bv.select(offset+1);
    for (uint64_t i = start; i < end; ++i) {
        seq += (char)seq_wt.at(i);
    }
    return seq;
}
    
/// Loop over all the handles to next/previous (right/left) nodes. Passes
/// them to a callback which returns false to stop iterating and true to
/// continue. Returns true if we finished and false if we stopped early.
bool graph_t::follow_edges(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const {
    bool result = true;
    bool offset = handle_helper::unpack_number(handle);
    bool is_rev = handle_helper::unpack_bit(handle);
    if (!go_left && !is_rev || go_left && is_rev) {
        uint64_t edges_begin = edge_fwd_wt.select(offset, 0);
        for (uint64_t i = edges_begin+1; ; ++i) {
            id_t id = edge_fwd_wt.at(i);
            if (!id) break; // end of record
            bool inv = edge_fwd_inv_bv.at(i);
            handle_t handle = handle_helper::pack(id, (inv ? !is_rev : is_rev));
            result &= iteratee(handle);
            if (!result) break;
        }
    } else {
        assert(go_left && !is_rev || !go_left && is_rev);
        uint64_t edges_begin = edge_rev_wt.select(offset, 0);
        for (uint64_t i = edges_begin+1; ; ++i) {
            id_t id = edge_rev_wt.at(i);
            if (!id) break; // end of record
            bool inv = edge_rev_inv_bv.at(i);
            handle_t handle = handle_helper::pack(id, (inv ? !is_rev : is_rev));
            result &= iteratee(handle);
            if (!result) break;
        }
    }
    return result;
}
    
/// Loop over all the nodes in the graph in their local forward
/// orientations, in their internal stored order. Stop if the iteratee
/// returns false. Can be told to run in parallel, in which case stopping
/// after a false return value is on a best-effort basis and iteration
/// order is not defined.
void graph_t::for_each_handle(const std::function<bool(const handle_t&)>& iteratee, bool parallel) const {
    if (parallel) {
        volatile bool flag=false;
#pragma omp parallel for
        for (uint64_t i = 0; i < graph_id_wt.size(); ++i) {
            if (flag) continue;
            bool result = iteratee(handle_helper::pack(i,false));
#pragma omp atomic
            flag &= result;
        }
    } else {
        for (uint64_t i = 0; i < graph_id_wt.size(); ++i) {
            if (!iteratee(handle_helper::pack(i,false))) break;
        }
    }
}
    
/// Return the number of nodes in the graph
/// TODO: can't be node_count because XG has a field named node_count.
size_t graph_t::node_size(void) const {
    graph_id_wt.size() - graph_id_wt.rank(graph_id_wt.size(), 0);
}
    
/// Return the smallest ID in the graph, or some smaller number if the
/// smallest ID is unavailable. Return value is unspecified if the graph is empty.
id_t graph_t::min_node_id(void) const {
    return _min_node_id;
}
    
/// Return the largest ID in the graph, or some larger number if the
/// largest ID is unavailable. Return value is unspecified if the graph is empty.
id_t graph_t::max_node_id(void) const {
    return _max_node_id;
}
    
////////////////////////////////////////////////////////////////////////////
// Interface that needs to be using'd
////////////////////////////////////////////////////////////////////////////
    
/// Loop over all the handles to next/previous (right/left) nodes. Works
/// with a callback that just takes all the handles and returns void.
/// Has to be a template because we can't overload on the types of std::function arguments otherwise.
/// MUST be pulled into implementing classes with `using` in order to work!
/*
template <typename T>
auto graph_t::follow_edges(const handle_t& handle, bool go_left, T&& iteratee) const
    -> typename std::enable_if<std::is_void<decltype(iteratee(get_handle(0, false)))>::value>::type {
    // Implementation only for void-returning iteratees
    // We ought to just overload on the std::function but that's not allowed until C++14.
    // See <https://stackoverflow.com/q/13811180>
        
    // We also can't use result_of<T(handle_t)>::type to sniff the return
    // type out because that ::type would not exist (since that's what you
    // get for a void apparently?) and we couldn't check if it's bool or
    // void.
        
    // So we do this nonsense thing with a trailing return type (to get the
    // actual arg into scope) and a decltype (which is allowed to resolve to
    // void) and is_void (which is allowed to take void) and a fake
    // get_handle call (which is the shortest handle_t-typed expression I
    // could think of).
        
    // Make a wrapper that puts a bool return type on.
    std::function<bool(const handle_t&)> lambda = [&](const handle_t& found) {
        iteratee(found);
        return true;
    };
        
    // Use that
    follow_edges(handle, go_left, lambda);
        
    // During development I managed to get earlier versions of this template to build infinitely recursive functions.
    static_assert(!std::is_void<decltype(lambda(get_handle(0, false)))>::value, "can't take our own lambda");
}
*/
    
/// Loop over all the nodes in the graph in their local forward
/// orientations, in their internal stored order. Works with void-returning iteratees.
/// MUST be pulled into implementing classes with `using` in order to work!
/*
template <typename T>
auto graph_t::for_each_handle(T&& iteratee, bool parallel) const
    -> typename std::enable_if<std::is_void<decltype(iteratee(get_handle(0, false)))>::value>::type {
    // Make a wrapper that puts a bool return type on.
    std::function<bool(const handle_t&)> lambda = [&](const handle_t& found) {
        iteratee(found);
        return true;
    };
        
    // Use that
    for_each_handle(lambda, parallel);
}
*/

/*
void graph_t::for_each_edge(const std::function<bool(const edge_t&)>& iteratee, bool parallel) {
    for_each_handle([&](const handle_t& handle){
            bool keep_going = true;
            // filter to edges where this node is lower ID or any rightward self-loops
            follow_edges(handle, false, [&](const handle_t& next) {
                    if (get_id(handle) <= get_id(next)) {
                        keep_going = iteratee(edge_handle(handle, next));
                    }
                    return keep_going;
                });
            if (keep_going) {
                // filter to edges where this node is lower ID or leftward reversing
                // self-loop
                follow_edges(handle, true, [&](const handle_t& prev) {
                        if (get_id(handle) < get_id(prev) ||
                            (get_id(handle) == get_id(prev) && !get_is_reverse(prev))) {
                            keep_going = iteratee(edge_handle(prev, handle));
                        }
                        return keep_going;
                    });
            }
        }, parallel);
}
*/
    
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
size_t graph_t::get_degree(const handle_t& handle, bool go_left) const {
    uint64_t offset = handle_helper::unpack_number(handle);
    bool is_rev = handle_helper::unpack_bit(handle);
    if (!go_left && !is_rev || go_left && is_rev) {
        return edge_fwd_wt.select(offset+1, 0) - edge_fwd_wt.select(offset, 0);
    } else {
        return edge_rev_wt.select(offset+1, 0) - edge_rev_wt.select(offset, 0);
    }
}
    
////////////////////////////////////////////////////////////////////////////
// Concrete utility methods
////////////////////////////////////////////////////////////////////////////
    
/// Get a Protobuf Visit from a handle.
//Visit to_visit(const handle_t& handle) const;
    
/// Get the locally forward version of a handle
handle_t graph_t::forward(const handle_t& handle) const {
    handle_t handle_fwd = handle;
    if (handle_helper::unpack_bit(handle)) handle_helper::toggle_bit(handle_fwd);
    return handle_fwd;
}

/*
/// A pair of handles can be used as an edge. When so used, the handles have a
/// canonical order and orientation.
edge_t graph_t::edge_handle(const handle_t& left, const handle_t& right) const {
    make_pair(left, right);
}
    
/// Such a pair can be viewed from either inward end handle and produce the
/// outward handle you would arrive at.
handle_t graph_t::traverse_edge_handle(const edge_t& edge, const handle_t& left) const {
}
*/
    
/**
 * This is the interface for a handle graph that stores embedded paths.
 */
    
////////////////////////////////////////////////////////////////////////////
// Path handle interface that needs to be implemented
////////////////////////////////////////////////////////////////////////////
    
/// Determine if a path name exists and is legal to get a path handle for.
bool graph_t::has_path(const std::string& path_name) const {
    std::string query_s = "$" + path_name + "$";
    std::vector<uint64_t> query_v(query_s.begin(), query_s.end());
    return path_name_fmi.locate(query_v).size() > 0;
}
    
/// Look up the path handle for the given path name.
/// The path with that name must exist.
path_handle_t graph_t::get_path_handle(const std::string& path_name) const {
    std::string query_s = "$" + path_name + "$";
    std::vector<uint64_t> query_v(query_s.begin(), query_s.end());
    std::vector<uint64_t> occs = path_name_fmi.locate(query_v);
    assert(occs.size() == 1);
    uint64_t offset = occs.front();
    return as_path_handle(path_name_bv.rank1(offset));
}
    
/// Look up the name of a path from a handle to it
std::string graph_t::get_path_name(const path_handle_t& path_handle) const {
    return paths.at(as_integer(path_handle)).name;
}
    
/// Returns the number of node occurrences in the path
size_t graph_t::get_occurrence_count(const path_handle_t& path_handle) const {
    return paths.at(as_integer(path_handle)).occurrence_count();
}

/// Returns the number of paths stored in the graph
size_t graph_t::get_path_count(void) const {
    return _path_count;
}
    
/// Execute a function on each path in the graph
// TODO: allow stopping early?
void graph_t::for_each_path_handle(const std::function<void(const path_handle_t&)>& iteratee) const {
    for (uint64_t i = 0; i < paths.size(); ++i) {
        if (paths.at(i).occurrence_count()) {
            iteratee(as_path_handle(i));
        }
    }
}
    
/// Get a node handle (node ID and orientation) from a handle to an occurrence on a path
handle_t graph_t::get_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    //path_handle_t = as_path_handle(occ_handle[0]);
    auto& path = paths.at(occ_handle[0]);
    // get the step
    step_t step = path.get_occurrence(occ_handle[1]);
    // compute the handle
    handle_t handle = as_handle(boundary_bv.rank1(step.start));
    // check if it's fully live
    size_t length = get_length(handle);
    for (uint64_t i = as_integer(handle); i < as_integer(handle)+length; ++i) {
        if (dead_wt.at(i) != 0) {
            return as_handle(0);
        }
    }
    // if it is, we return it
    return handle;
}

/// Get a handle to the first occurrence in a path.
/// The path MUST be nonempty.
occurrence_handle_t graph_t::get_first_occurrence(const path_handle_t& path_handle) const {
    auto& path = paths.at(as_integer(path_handle));
    assert(path.occurrence_count());
    occurrence_handle_t result;
    int64_t* r_ints = as_integers(result);
    r_ints[0] = as_integer(path_handle);
    r_ints[1] = 0;
    return result;
}
    
/// Get a handle to the last occurrence in a path
/// The path MUST be nonempty.
occurrence_handle_t graph_t::get_last_occurrence(const path_handle_t& path_handle) const {
    auto& path = paths.at(as_integer(path_handle));
    assert(path.occurrence_count());
    occurrence_handle_t result;
    int64_t* r_ints = as_integers(result);
    r_ints[0] = as_integer(path_handle);
    r_ints[1] = path.occurrence_count()-1;
    return result;
}
    
/// Returns true if the occurrence is not the last occurence on the path, else false
bool graph_t::has_next_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    auto& path = paths.at(occ_handle[0]);
    return occ_handle[1] < path.occurrence_count()-1;
}
    
/// Returns true if the occurrence is not the first occurence on the path, else false
bool graph_t::has_previous_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    auto& path = paths.at(occ_handle[0]);
    assert(path.occurrence_count());
    return occ_handle[1] > 0;
}

/// Returns a handle to the next occurrence on the path
occurrence_handle_t graph_t::get_next_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    auto& path = paths.at(occ_handle[0]);
    assert(path.occurrence_count());
    uint64_t rank = occ_handle[1];
    occurrence_handle_t next_occ_handle = occurrence_handle;
    int64_t* next_i = as_integers(next_occ_handle);
    next_i[1] = rank+1;
    return next_occ_handle;
}

/// Returns a handle to the previous occurrence on the path
occurrence_handle_t graph_t::get_previous_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    auto& path = paths.at(occ_handle[0]);
    assert(path.occurrence_count());
    uint64_t rank = occ_handle[1];
    occurrence_handle_t prev_occ_handle = occurrence_handle;
    int64_t* prev_i = as_integers(prev_occ_handle);
    prev_i[1] = rank-1;
    return prev_occ_handle;
}
    
/// Returns a handle to the path that an occurrence is on
path_handle_t graph_t::get_path_handle_of_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    return as_path_handle(occ_handle[0]);
}
    
/// Returns the 0-based ordinal rank of a occurrence on a path
size_t graph_t::get_ordinal_rank_of_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    return occ_handle[1];
}

////////////////////////////////////////////////////////////////////////////
// Additional optional interface with a default implementation
////////////////////////////////////////////////////////////////////////////

/// Returns true if the given path is empty, and false otherwise
bool graph_t::is_empty(const path_handle_t& path_handle) const {
    auto& path = paths.at(as_integer(path_handle));
    return path.occurrence_count() == 0;
}

////////////////////////////////////////////////////////////////////////////
// Concrete utility methods
////////////////////////////////////////////////////////////////////////////

/// Loop over all the occurrences along a path, from first through last
void graph_t::for_each_occurrence_in_path(const path_handle_t& path, const std::function<void(const occurrence_handle_t&)>& iteratee) const {
    if (is_empty(path)) return;
    occurrence_handle_t occ = get_first_occurrence(path);
    iteratee(occ); // run the first time
    while (has_next_occurrence(occ)) {
        occ = get_next_occurrence(occ);
        iteratee(occ);
    }
}

/**
 * This is the interface for a handle graph that supports modification.
 */
/*
 * Note: All operations may invalidate path handles and occurrence handles.
 */
    
/// Create a new node with the given sequence and return the handle.
handle_t graph_t::create_handle(const std::string& sequence) {
}

/// Create a new node with the given id and sequence, then return the handle.
handle_t graph_t::create_handle(const std::string& sequence, const id_t& id) {
}
    
/// Remove the node belonging to the given handle and all of its edges.
/// Does not update any stored paths.
/// Invalidates the destroyed handle.
/// May be called during serial for_each_handle iteration **ONLY** on the node being iterated.
/// May **NOT** be called during parallel for_each_handle iteration.
/// May **NOT** be called on the node from which edges are being followed during follow_edges.
void graph_t::destroy_handle(const handle_t& handle) {
}
    
/// Create an edge connecting the given handles in the given order and orientations.
/// Ignores existing edges.
void graph_t::create_edge(const handle_t& left, const handle_t& right) {
}
    
/// Remove the edge connecting the given handles in the given order and orientations.
/// Ignores nonexistent edges.
/// Does not update any stored paths.
void graph_t::destroy_edge(const handle_t& left, const handle_t& right) {
}
        
/// Remove all nodes and edges. Does not update any stored paths.
void graph_t::clear(void) {
}
    
/// Swap the nodes corresponding to the given handles, in the ordering used
/// by for_each_handle when looping over the graph. Other handles to the
/// nodes being swapped must not be invalidated. If a swap is made while
/// for_each_handle is running, it affects the order of the handles
/// traversed during the current traversal (so swapping an already seen
/// handle to a later handle's position will make the seen handle be visited
/// again and the later handle not be visited at all).
void graph_t::swap_handles(const handle_t& a, const handle_t& b) {
}
    
/// Alter the node that the given handle corresponds to so the orientation
/// indicated by the handle becomes the node's local forward orientation.
/// Rewrites all edges pointing to the node and the node's sequence to
/// reflect this. Invalidates all handles to the node (including the one
/// passed). Returns a new, valid handle to the node in its new forward
/// orientation. Note that it is possible for the node's ID to change.
/// Does not update any stored paths. May change the ordering of the underlying
/// graph.
handle_t graph_t::apply_orientation(const handle_t& handle) {
}
    
/// Split a handle's underlying node at the given offsets in the handle's
/// orientation. Returns all of the handles to the parts. Other handles to
/// the node being split may be invalidated. The split pieces stay in the
/// same local forward orientation as the original node, but the returned
/// handles come in the order and orientation appropriate for the handle
/// passed in.
/// Updates stored paths.
std::vector<handle_t> graph_t::divide_handle(const handle_t& handle, const std::vector<size_t>& offsets) {
}
    

/**
 * This is the interface for a handle graph with embedded paths where the paths can be modified.
 * Note that if the *graph* can also be modified, the implementation will also
 * need to inherit from MutableHandleGraph, via the combination
 * MutablePathMutableHandleGraph interface.
 * TODO: This is a very limited interface at the moment. It will probably need to be extended.
 */
    
/**
 * Destroy the given path. Invalidates handles to the path and its node occurrences.
 */
void graph_t::destroy_path(const path_handle_t& path) {
}

/**
 * Create a path with the given name. The caller must ensure that no path
 * with the given name exists already, or the behavior is undefined.
 * Returns a handle to the created empty path. Handles to other paths must
 * remain valid.
 */
path_handle_t graph_t::create_path_handle(const std::string& name) {
}
    
/**
 * Append a visit to a node to the given path. Returns a handle to the new
 * final occurrence on the path which is appended. Handles to prior
 * occurrences on the path, and to other paths, must remain valid.
 */
occurrence_handle_t graph_t::append_occurrence(const path_handle_t& path, const handle_t& to_append) {
}


}
