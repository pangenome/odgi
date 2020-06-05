//
//  odgi.cpp
//  

#include "odgi.hpp"

namespace odgi {

void ODGI::serialize_members(std::ostream& out) const {
    serialize_and_measure(out);
}

void ODGI::deserialize_members(std::istream& in) {
    load(in);
}

uint32_t ODGI::get_magic_number() const {
    return 2880272692ul;
}

/// Method to check if a node exists by ID
bool ODGI::has_node(nid_t node_id) const {
    if (graph_id_hidden_set.count(node_id)) {
        // Hide this node.
        return false;
    }
    uint64_t rank = id_to_rank(node_id);
    return (rank >= node_v.size() ? false : !deleted_node_bv.at(rank));
}

/// Look up the handle for the node with the given ID in the given orientation
handle_t ODGI::get_handle(const nid_t& node_id, bool is_reverse) const {
    return number_bool_packing::pack(id_to_rank(node_id), is_reverse);
}

/// Get the ID from a handle
nid_t ODGI::get_id(const handle_t& handle) const {
    return rank_to_id(number_bool_packing::unpack_number(handle));
}

/// get the backing node rank for a given node id
uint64_t ODGI::id_to_rank(const nid_t& node_id) const {
    return node_id - _id_increment;
}

/// get the node id for the node with the given backing rank
nid_t ODGI::rank_to_id(const uint64_t& rank) const {
    return rank + _id_increment;
}

/// set the id increment, used when the graph starts at a high id to reduce loading costs
void ODGI::set_id_increment(const nid_t& min_id) {
    if (min_id < 1) {
        // Clamp at the low end to 1, because rank 0 needs at least ID 1.
        set_id_increment(1);
        return;
    }
    
    if ((get_node_count() + graph_id_hidden_set.size()) > 0) {
        increment_node_ids(min_id - _id_increment);
    }
    _id_increment = min_id;
}

/// Add the given value to all node IDs
void ODGI::increment_node_ids(nid_t increment) {
    if (deleted_node_bv.size() == _deleted_node_count) {
        // Nothing is in this graph! Nothing to do!
        return;
    }

    if (increment > 0) {
        // Shift right
        
        // Make sure node and deletion vectors are big enough
        size_t old_size = node_v.size();
        size_t new_size = old_size + increment;
        node_v.resize(new_size);
        while (deleted_node_bv.size() < new_size + 1) {
            // Put a bunch of new 0s before the final delimiting 1.
            deleted_node_bv.insert(old_size, 0);
        }
        // Shift all the nodes and deletion flags up
        for (size_t i = new_size - 1; i >= increment; i--) {
            node_v[i] = node_v[i - increment];
            deleted_node_bv[i] = deleted_node_bv[i - increment];
        }
        // Clear out the low end (everyone is deleted)
        for (size_t i = 0; i < increment; i++) {
            deleted_node_bv[i] = 1;
        }
        _deleted_node_count += increment;
    } else if (increment < 0) {
        // Shift left
        
        if (-increment >= node_v.size()) {
            throw std::runtime_error("Cannot shift IDs left by " + std::to_string(-increment) +
                " because graph nodes would be destroyed.");
        }
        
        for (size_t i = 0; i < -increment; i++) {
            // Make sure all the slots we are clobbering are empty.
            if (!deleted_node_bv[i]) {
                throw std::runtime_error("Cannot shift IDs left by " + std::to_string(-increment) +
                    " because graph nodes would be destroyed.");
            }
        }
        
        // Shift everything down
        for (size_t i = 0; i < node_v.size() + increment; i++) {
            node_v[i] = node_v[i - increment];
            deleted_node_bv[i] = deleted_node_bv[i - increment];
        }
        // Clear out the high end (everyone is deleted).
        // No need to touch the final delimiting 1.
        for (size_t i = deleted_node_bv.size() + increment; i < deleted_node_bv.size() - 1; i++) {
            deleted_node_bv[i] = 1;
        }
        
        // Total number of deleted nodes is conserved.
    } else {
        // No shift. Nothing to do!
        return;
    }
        
    // Replace the hidden set
    decltype(graph_id_hidden_set) new_hidden_set;
    for (auto& item : graph_id_hidden_set) {
        new_hidden_set.insert(item + increment);
    }
    graph_id_hidden_set = std::move(new_hidden_set);
    
    // Adjust the min and max node ranks
    _max_node_rank += increment;
    _min_node_rank += increment;

    for (auto& path : path_metadata_v) {
        if (path.length > 0) {
            // Increment all the ranks stored in the normal handles inside the step
            // handles for path begin and end.
            set_handle_of_step(path.first, add_to_number(get_handle_of_step(path.first), increment));
            set_handle_of_step(path.last, add_to_number(get_handle_of_step(path.last), increment));
        }
    }
}

/// Reassign all node IDs as specified by the old->new mapping function.
void ODGI::reassign_node_ids(const std::function<nid_t(const nid_t&)>& get_new_id) {
    ODGI reassigned;
    // nodes
    for_each_handle(
        [&](const handle_t& handle) {
            reassigned.create_handle(get_sequence(handle), get_new_id(get_id(handle)));
        });
    // edges
    for_each_edge(
        [&](const edge_t& edge) {
            reassigned.create_edge(reassigned.get_handle(get_new_id(get_id(edge.first)), get_is_reverse(edge.first)),
                                   reassigned.get_handle(get_new_id(get_id(edge.second)), get_is_reverse(edge.second)));
        });
    // paths
    for_each_path_handle(
        [&](const path_handle_t& old_path) {
            path_handle_t new_path = reassigned.create_path_handle(get_path_name(old_path));
            for_each_step_in_path(
                old_path,
                [&](const step_handle_t& step) {
                    handle_t old_handle = get_handle_of_step(step);
                    handle_t new_handle = reassigned.get_handle(get_new_id(get_id(old_handle)),
                                                                get_is_reverse(old_handle));
                    reassigned.append_step(new_path, new_handle);
                });
            if (get_is_circular(old_path)) {
                reassigned.set_circularity(new_path, true);
            }
        });
    *this = reassigned;
}
    
/// Get the orientation of a handle
bool ODGI::get_is_reverse(const handle_t& handle) const {
    return number_bool_packing::unpack_bit(handle);
}
    
/// Invert the orientation of a handle (potentially without getting its ID)
handle_t ODGI::flip(const handle_t& handle) const {
    return number_bool_packing::toggle_bit(handle);
}
    
/// Get the length of a node
size_t ODGI::get_length(const handle_t& handle) const {
    return node_v.at(number_bool_packing::unpack_number(handle)).sequence_size();
}

/// Get the sequence of a node, presented in the handle's local forward orientation.
std::string ODGI::get_sequence(const handle_t& handle) const {
    auto& seq = node_v.at(number_bool_packing::unpack_number(handle)).sequence();
    return (get_is_reverse(handle) ? reverse_complement(seq) : seq);
}
    
/// Loop over all the handles to next/previous (right/left) nodes. Passes
/// them to a callback which returns false to stop iterating and true to
/// continue. Returns true if we finished and false if we stopped early.
bool ODGI::follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const {
    const node_t& node = node_v.at(number_bool_packing::unpack_number(handle));
    bool is_rev = get_is_reverse(handle);
    nid_t node_id = get_id(handle);
    const std::vector<uint64_t> node_edges = node.edges();
    if (node_edges.size() == 0) return true;
    for (uint64_t i = 0; i < node_edges.size(); i+=2) {
        // unpack the edge
        uint64_t other_id = edge_delta_to_id(node_id, node_edges.at(i));
        uint8_t packed_edge = node_edges.at(i+1);
        bool on_rev = edge_helper::unpack_on_rev(packed_edge);
        bool other_rev = edge_helper::unpack_other_rev(packed_edge);
        bool to_curr = edge_helper::unpack_to_curr(packed_edge);
        if (other_id == node_id && on_rev == other_rev) {
            // non-inverting self loop
            // we can go either direction
            to_curr = go_left;
            other_rev = is_rev;
        } else if (is_rev != on_rev) {
            other_rev ^= 1;
            to_curr ^= 1;
        }
        if (!go_left && !to_curr
            || go_left && to_curr) {
            if (!iteratee(get_handle(other_id, other_rev))) {
                return false;
            }
        }
    }
    return true;
}

    
/// Loop over all the nodes in the graph in their local forward
/// orientations, in their internal stored order. Stop if the iteratee
/// returns false. Can be told to run in parallel, in which case stopping
/// after a false return value is on a best-effort basis and iteration
/// order is not defined.
bool ODGI::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel) const {
    if (parallel) {
        volatile bool flag=true;
#pragma omp parallel for
        for (uint64_t i = 0; i < node_v.size(); ++i) {
            if (deleted_node_bv.at(i) == 1) continue;
            if (!flag) continue;
            handle_t packed = number_bool_packing::pack(i,false);
            if (graph_id_hidden_set.count(get_id(packed))) continue;
            bool result = iteratee(packed);
#pragma omp atomic
            flag &= result;
        }
        return flag;
    } else {
        for (uint64_t i = 0; i < node_v.size(); ++i) {
            if (deleted_node_bv.at(i) == 1) continue;
            handle_t packed = number_bool_packing::pack(i,false);
            if (graph_id_hidden_set.count(get_id(packed))) continue;
            if (!iteratee(packed)) return false;
        }
        return true;
    }
}

/// Return the number of nodes in the graph.
/// Doesn't count deleted or hidden nodes.
/// TODO: can't be node_count because XG has a field named node_count.
size_t ODGI::get_node_count(void) const {
    return node_v.size() - _deleted_node_count - graph_id_hidden_set.size();
}
    
/// Return the smallest ID in the graph, or some smaller number if the
/// smallest ID is unavailable. Return value is unspecified if the graph is empty.
nid_t ODGI::min_node_id(void) const {
    // May overflow for empty graphs.
    return rank_to_id(_min_node_rank);
}
    
/// Return the largest ID in the graph, or some larger number if the
/// largest ID is unavailable. Return value is unspecified if the graph is empty.
nid_t ODGI::max_node_id(void) const {
    return rank_to_id(_max_node_rank);
}
    
////////////////////////////////////////////////////////////////////////////
// Additional optional interface with a default implementation
////////////////////////////////////////////////////////////////////////////
    
/// Get the number of edges on the right (go_left = false) or left (go_left
/// = true) side of the given handle. The default implementation is O(n) in
/// the number of edges returned, but graph implementations that track this
/// information more efficiently can override this method.
size_t ODGI::get_degree(const handle_t& handle, bool go_left) const {
    size_t degree = 0;
    follow_edges(handle, go_left, [&degree](const handle_t& h) { ++degree; });
    return degree;
}


/// Get the locally forward version of a handle
handle_t ODGI::forward(const handle_t& handle) const {
    return HandleGraph::forward(handle);
}

/// A pair of handles can be used as an edge. When so used, the handles have a
/// canonical order and orientation.
edge_t ODGI::edge_handle(const handle_t& left, const handle_t& right) const {
    return HandleGraph::edge_handle(left, right);
}

/// View the given edge handle from either inward end handle and produce the
/// outward handle you would arrive at.
handle_t ODGI::traverse_edge_handle(const edge_t& edge, const handle_t& left) const {
    return HandleGraph::traverse_edge_handle(edge, left);
}
    
/**
 * This is the interface for a handle graph that stores embedded paths.
 */
    
////////////////////////////////////////////////////////////////////////////
// Path handle interface that needs to be implemented
////////////////////////////////////////////////////////////////////////////
    
/// Determine if a path name exists and is legal to get a path handle for.
bool ODGI::has_path(const std::string& path_name) const {
    auto f = path_name_map.find(path_name);
    if (f == path_name_map.end()) return false;
    else return true;
}
    
/// Look up the path handle for the given path name.
/// The path with that name must exist.
path_handle_t ODGI::get_path_handle(const std::string& path_name) const {
    auto f = path_name_map.find(path_name);
    assert(f != path_name_map.end());
    return f->second;
}

/// Look up the name of a path from a handle to it
std::string ODGI::get_path_name(const path_handle_t& path_handle) const {
    return find_metadata(path_handle).name;
}
    
/// Returns the number of node steps in the path
size_t ODGI::get_step_count(const path_handle_t& path_handle) const {
    return find_metadata(path_handle).length;
}

/// Returns the number of paths stored in the graph
size_t ODGI::get_path_count(void) const {
    return _path_count;
}
    
/// Execute a function on each path in the graph
bool ODGI::for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const {
    bool flag = true;
    for (uint64_t i = 0; i < _path_rank_next && flag; ++i) {
        path_handle_t path = rank_to_path(i);
        if (get_step_count(path) > 0) {
            flag &= iteratee(path);
        }
    }
    return flag;
}

bool ODGI::for_each_step_on_handle_impl(const handle_t& handle, const std::function<bool(const step_handle_t&)>& iteratee) const {
    uint64_t handle_n = number_bool_packing::unpack_number(handle);
    const node_t& node = node_v.at(handle_n);
    const std::vector<node_t::step_t> steps = node.get_path_steps();
    bool is_rev = get_is_reverse(handle);
    bool flag = true;
    uint64_t i = 0;
    for (auto& step : steps) {
        step_handle_t step_handle;
        as_integers(step_handle)[0] = as_integer(number_bool_packing::pack(handle_n, step.is_rev()));
        as_integers(step_handle)[1] = i++;
        flag &= iteratee(step_handle);
    }
    return flag;
}

/// Returns a vector of all steps of a node on paths. Optionally restricts to
/// steps that match the handle in orientation.
std::vector<step_handle_t> ODGI::steps_of_handle(const handle_t& handle,
                                                    bool match_orientation) const {
    std::vector<step_handle_t> res;
    for_each_step_on_handle(handle, [&](const step_handle_t& step) {
            handle_t h = get_handle_of_step(step);
            if (!match_orientation || get_is_reverse(h) == get_is_reverse(handle)) {
                res.push_back(step);
            }
        });
    return res;
}

size_t ODGI::get_step_count(const handle_t& handle) const {
    const node_t& node = node_v.at(number_bool_packing::unpack_number(handle));
    return node.path_count();
}

/// Get a node handle (node ID and orientation) from a handle to an step on a path
handle_t ODGI::get_handle_of_step(const step_handle_t& step_handle) const {
    assert(!is_path_front_end(step_handle));
    assert(!is_path_end(step_handle));
    return as_handle(as_integers(step_handle)[0]);
}

/// Get a handle to the first step in a path.
/// The path MUST be nonempty.
step_handle_t ODGI::path_begin(const path_handle_t& path_handle) const {
    return find_metadata(path_handle).first;
}
    
/// Get a handle to the last step in a path
/// The path MUST be nonempty.
step_handle_t ODGI::path_back(const path_handle_t& path_handle) const {
    return find_metadata(path_handle).last;
}

/// Get the reverse end iterator
step_handle_t ODGI::path_front_end(const path_handle_t& path_handle) const {
    step_handle_t step;
    as_integers(step)[0] = path_to_rank(path_handle);
    as_integers(step)[1] = std::numeric_limits<uint64_t>::max() - 1;
    return step;
}

/// Get the forward end iterator
step_handle_t ODGI::path_end(const path_handle_t& path_handle) const {
    step_handle_t step;
    as_integers(step)[0] = path_to_rank(path_handle);
    as_integers(step)[1] = std::numeric_limits<uint64_t>::max();
    return step;
}

/// is the path circular
bool ODGI::get_is_circular(const path_handle_t& path_handle) const {
    return find_metadata(path_handle).is_circular;
}

/// set the circular flag for the path
void ODGI::set_circularity(const path_handle_t& path_handle, bool circular) {
    find_metadata(path_handle).is_circular = circular;
}
    
/// Returns true if the step has a next step on the path that doesn't wrap around, else false
bool ODGI::has_linear_next_step(const step_handle_t& step_handle) const {
    if (is_path_end(step_handle)) {
        // This is a special end handle with nowhere to go next
        return false;
    }
    if (is_path_front_end(step_handle)) {
        // This is the front-end handle; we only have a successor if the path is nonempty
        return get_step_count(get_path_handle_of_step(step_handle)) > 0; 
    }
    // Otherwise we have a node, so look it up.
    const node_t& node = node_v.at(number_bool_packing::unpack_number(get_handle_of_step(step_handle)));
    return node.get_path_step(as_integers(step_handle)[1]).next_id() != path_end_marker;
}
    
/// Returns true if the step has a previous step on the path that doesn't wrap around, else false
bool ODGI::has_linear_previous_step(const step_handle_t& step_handle) const {
    if (is_path_front_end(step_handle)) {
        // This is a special front-end handle with nowhere to come from
        return false;
    }
    if (is_path_end(step_handle)) {
        // This is the end handle; we only have a predecessor if the path is nonempty
        return get_step_count(get_path_handle_of_step(step_handle)) > 0; 
    }
    // Otherwise we have a node, so look it up.
    const node_t& node = node_v.at(number_bool_packing::unpack_number(get_handle_of_step(step_handle)));
    return node.get_path_step(as_integers(step_handle)[1]).prev_id() != path_begin_marker;
}
    
/// Returns true if the step has a next step on the path, else false
bool ODGI::has_next_step(const step_handle_t& step_handle) const {
    return has_linear_next_step(step_handle) || get_is_circular(get_path_handle_of_step(step_handle));
}
    
/// Returns true if the step has a previous step on the path, else false
bool ODGI::has_previous_step(const step_handle_t& step_handle) const {
    return has_linear_previous_step(step_handle) || get_is_circular(get_path_handle_of_step(step_handle));
}

bool ODGI::is_path_front_end(const step_handle_t& step_handle) const {
    return as_integers(step_handle)[1] == std::numeric_limits<uint64_t>::max() - 1;
}

bool ODGI::is_path_end(const step_handle_t& step_handle) const {
    return as_integers(step_handle)[1] == std::numeric_limits<uint64_t>::max();
}

/// Returns a handle to the next step on the path
/// Returns the forward end iterator if none exists
step_handle_t ODGI::get_next_step(const step_handle_t& step_handle) const {
    // check if we have a magic path iterator step handle
    if (is_path_front_end(step_handle)) {
        // After front end is the beginning of the path. Path is stored as the step's first number.
        return path_begin(rank_to_path(as_integers(step_handle)[0]));
    } else if (is_path_end(step_handle)) {
        // Nothing really comes after this. Caller probably shouldn't be calling us. Just stay here.
        return step_handle;
    }
    // Otherwise we're actually on a node.
    handle_t curr_handle = get_handle_of_step(step_handle);
    nid_t curr_id = get_id(curr_handle);
    const node_t& node = node_v.at(number_bool_packing::unpack_number(curr_handle));
    auto& step = node.get_path_step(as_integers(step_handle)[1]);
    if (step.next_id() == path_end_marker) {
        if (find_metadata(rank_to_path(step.path_id())).is_circular) {
            // End of a circular path wraps back to start
            return find_metadata(rank_to_path(step.path_id())).first;
        } else {
            return path_end(get_path_handle_of_step(step_handle));
        }
    }
    nid_t next_id = edge_delta_to_id(curr_id, step.next_id() - (path_end_marker + 1));
    handle_t next_handle = get_handle(next_id);
    bool next_rev = node_v.at(number_bool_packing::unpack_number(next_handle)).get_path_step(step.next_rank()).is_rev();
    step_handle_t next_step;
    as_integers(next_step)[0] = as_integer(get_handle(next_id, next_rev));
    as_integers(next_step)[1] = step.next_rank();
    return next_step;
}

/// Returns a handle to the previous step on the path
step_handle_t ODGI::get_previous_step(const step_handle_t& step_handle) const {
    // check if we have a magic path iterator step handle
    if (is_path_front_end(step_handle)) {
        // Nothing really comes before this. Caller probably shouldn't be calling us. Just stay here.
        return step_handle;
    } else if (is_path_end(step_handle)) {
        // Before end is the last thing on the path. Path is stored as the step's first number.
        return path_back(rank_to_path(as_integers(step_handle)[0]));
    }
    // Otherwise we're actually on a node.
    handle_t curr_handle = get_handle_of_step(step_handle);
    //handle_t curr_handle = get_handle_of_step(step_handle);
    nid_t curr_id = get_id(curr_handle);
    const node_t& node = node_v.at(number_bool_packing::unpack_number(curr_handle));
    auto& step = node.get_path_step(as_integers(step_handle)[1]);
    if (step.prev_id() == path_begin_marker) {
        if (find_metadata(rank_to_path(step.path_id())).is_circular) {
            // End of a start path wraps back to end
            return find_metadata(rank_to_path(step.path_id())).last;
        } else {
            return path_front_end(get_path_handle_of_step(step_handle));
        }
    }
    nid_t prev_id = edge_delta_to_id(curr_id, step.prev_id() - (path_end_marker + 1));
    handle_t prev_handle = get_handle(prev_id);
    bool prev_rev = node_v.at(number_bool_packing::unpack_number(prev_handle)).get_path_step(step.prev_rank()).is_rev();
    step_handle_t prev_step;
    as_integers(prev_step)[0] = as_integer(get_handle(prev_id, prev_rev));
    as_integers(prev_step)[1] = step.prev_rank();
    return prev_step;
}

path_handle_t ODGI::get_path_handle_of_step(const step_handle_t& step_handle) const {
    if (is_path_front_end(step_handle) || is_path_end(step_handle)) {
        // Looks like a sentinel. Path is stored as the first number.
        return rank_to_path(as_integers(step_handle)[0]);
    }
    
    // Otherwise, this step actually visits a node, so we ask the node what path it is on.
    const node_t& node = node_v.at(number_bool_packing::unpack_number(get_handle_of_step(step_handle)));
    auto& step = node.get_path_step(as_integers(step_handle)[1]);
    return rank_to_path(step.path_id());
}
    
////////////////////////////////////////////////////////////////////////////
// Additional optional interface with a default implementation
////////////////////////////////////////////////////////////////////////////

/// Returns true if the given path is empty, and false otherwise
bool ODGI::is_empty(const path_handle_t& path_handle) const {
    return get_step_count(path_handle) == 0;
}

/**
 * This is the interface for a handle graph that supports modification.
 */
/*
 * Note: All operations may invalidate path handles and step handles.
 */

/// Loop over all the steps along a path, from first through last
void ODGI::for_each_step_in_path(const path_handle_t& path, const std::function<void(const step_handle_t&)>& iteratee) const {
    auto& p = find_metadata(path);
    if (is_empty(path)) return;
    bool is_circular = get_is_circular(path);
    step_handle_t begin_step = path_begin(path);
    step_handle_t step = begin_step; // copy
    bool keep_going = true;
    do {
        iteratee(step);
        // if it's circular, and there is a next step that isn't the same as the beginning, continue
        // if it's not circular, and there is a next step, continue
        if (is_circular) {
            step = get_next_step(step);
            if (step == begin_step) {
                keep_going = false;
            }
        } else {
            if (has_next_step(step)) {
                step = get_next_step(step);
            } else {
                keep_going = false;
            }
        }
    } while (keep_going);
}
    
/// Create a new node with the given sequence and return the handle.
handle_t ODGI::create_handle(const std::string& sequence) {
    // get first deleted node to recycle
    if (_deleted_node_count) {
        uint64_t replaced_rank = deleted_node_bv.select1(0);
        nid_t replaced_id = rank_to_id(replaced_rank);
        return create_handle(sequence, replaced_id);
    } else {
        return create_handle(sequence, rank_to_id(node_v.size()));
    }
}

handle_t ODGI::create_hidden_handle(const std::string& sequence) {
    // get node id as max+1
    handle_t handle = create_handle(sequence);
    nid_t id = get_id(handle);
    graph_id_hidden_set.insert(id);
    return handle;
}

/// Create a new node with the given id and sequence, then return the handle.
handle_t ODGI::create_handle(const std::string& sequence, const nid_t& id) {
    
    if (sequence.empty()) {
        throw std::runtime_error("error:[ODGI] tried to create an empty node with ID " + std::to_string(id));
    }
    
    if (id <= 0) {
        throw std::runtime_error("error:[ODGI] tried to create a node with non-positive ID " + std::to_string(id));
    }
    
    if (id < _id_increment) {
        // We would need to change the ID increment to take this node, and
        // invalidate existing handles, which create_handle can't do.
        throw std::runtime_error("ODGI graph cannot accept a new node with ID " +
            std::to_string(id) + " in a graph with ID increment " + 
            std::to_string(_id_increment));
    }
    
    // Give the node a rank, which is just how far above _id_increment it is.
    uint64_t handle_rank = id_to_rank(id);
    
    if (handle_rank >= node_v.size()) {
        uint64_t old_size = node_v.size();
        // realloc to have a last entry handle_rank
        node_v.resize(handle_rank + 1);
        _node_count = node_v.size();
        
        // mark empty nodes (treat them as deleted)
        uint64_t added = _node_count - old_size;
        for (uint64_t i = 0; i < added; ++i) {
            // insert before final delimiter
            deleted_node_bv.insert(old_size, 1);
            ++_deleted_node_count;
        }
    }
    // update min/max node internal ranks
    _max_node_rank = max(handle_rank, _max_node_rank);
    _min_node_rank = min(handle_rank, _min_node_rank);
    
    // add to node vector
    // set its values
    auto& node = node_v[handle_rank];
    node.set_sequence(sequence);
    // it's not deleted
    assert(deleted_node_bv[handle_rank]);
    deleted_node_bv[handle_rank] = 0;
    --_deleted_node_count;
    // return handle
    return number_bool_packing::pack(handle_rank, 0);
}
    
/// Remove the node belonging to the given handle and all of its edges.
/// If any paths visit it, it becomes a "hidden" node accessible only via the paths.
void ODGI::destroy_handle(const handle_t& handle) {
    handle_t fwd_handle = get_is_reverse(handle) ? flip(handle) : handle;
    uint64_t id = get_id(handle);
    // remove steps in edge lists
    // enumerate the edges
    std::vector<edge_t> edges_to_destroy;
    follow_edges(fwd_handle, false, [&edges_to_destroy,&fwd_handle,this](const handle_t& h) {
            edges_to_destroy.push_back(make_pair(fwd_handle, h)); });
    follow_edges(fwd_handle, true, [&edges_to_destroy,&fwd_handle,this](const handle_t& h) {
            edges_to_destroy.push_back(make_pair(h, fwd_handle)); });
    // and then remove them
    for (auto& edge : edges_to_destroy) {
        destroy_edge(edge);
    }
   
    // Work out if it has any paths on it
    auto& node = node_v[number_bool_packing::unpack_number(handle)];
    if (node.path_count() > 0) {
        // Paths visit here, so hide the handle.
        // We will get called again when the last one leaves.
        graph_id_hidden_set.insert(id);
    } else {
        // No paths visit here. Destroy the handle and make sure it isn't hidden if it was.
        
        // clear the node storage
        node.clear();
        // remove from the graph by hiding it (compaction later)
        deleted_node_bv[number_bool_packing::unpack_number(handle)] = 1;
        // and from the set of hidden nodes, if it's a member
        if (graph_id_hidden_set.count(id)) {
            graph_id_hidden_set.erase(id);
        }
        ++_deleted_node_count;
        // TODO: check if we should compact our deleted nodes storage
    }
    
    
}

/*
void ODGI::rebuild_id_handle_mapping(void) {
    // for each live node, record the id in a new vector
    uint64_t j = 0;
    for (uint64_t i = 0; i < node_v.size(); ++i) {
        auto id = node_v[i].id();
        if (id == 0) continue;
        graph_id_map[id] = j++;
    }
    for (uint64_t i = 0; i < node_v.size(); ++i) {
        auto id = node_v[i].id();
        if (id == 0) {
            deleted_node_bv.remove(i);
        }
    }
    node_v.erase(std::remove_if(node_v.begin(), node_v.end(),
                                [&](const node_t& node) { return node.id() == 0; }),
                 node_v.end());
    _deleted_node_count = 0;
}
*/
    
/// Create an edge connecting the given handles in the given order and orientations.
/// Ignores existing edges.
void ODGI::create_edge(const handle_t& left_h, const handle_t& right_h) {
    //if (has_edge(left, right)) return; // do nothing if edge exists
    /*
    std::cerr << "create_edge " << get_id(left_h) << ":" << get_is_reverse(left_h)
              << " -> "
              << get_id(right_h) << ":" << get_is_reverse(right_h) << std::endl;
    */
    if (has_edge(left_h, right_h)) return; // do nothing if edge exists

    uint64_t left_rank = number_bool_packing::unpack_number(left_h);
    uint64_t right_rank = number_bool_packing::unpack_number(right_h);
    uint64_t right_relative = edge_to_delta(right_h, left_h);
    uint64_t left_relative = edge_to_delta(left_h, right_h);
    
    // insert the edge for each side
    auto& left_node = node_v.at(left_rank);
    left_node.add_edge(left_relative,
                       edge_helper::pack(get_is_reverse(left_h),
                                         get_is_reverse(right_h),
                                         false));

    ++_edge_count;
    // only insert the second side if it's on a different node
    if (left_rank == right_rank) return;

    auto& right_node = node_v.at(right_rank);
    right_node.add_edge(right_relative,
                        edge_helper::pack(get_is_reverse(right_h),
                                          get_is_reverse(left_h),
                                          true));

}

uint64_t ODGI::edge_delta_to_id(uint64_t base, uint64_t delta) const {
    assert(delta != 0);
    if (delta == 1) {
        return base;
    } else if (delta % 2 == 0) {
        return base + delta/2;
    } else { //if (delta-1 % 2 == 0) {
        return base - (delta-1)/2;
    }
}

uint64_t ODGI::edge_to_delta(const handle_t& left, const handle_t& right) const {
    int64_t delta = get_id(right) - get_id(left);
    return (delta == 0 ? 1 : (delta > 0 ? 2*abs(delta) : 2*abs(delta)+1));
}

bool ODGI::has_edge(const handle_t& left, const handle_t& right) const {
    bool exists = false;
    follow_edges(left, false, [&right, &exists](const handle_t& next) {
            if (next == right) exists = true;
            return !exists;
        });
    return exists;
}

/// Remove the edge connecting the given handles in the given order and orientations.
/// Ignores nonexistent edges.
/// Does not update any stored paths.
void ODGI::destroy_edge(const handle_t& left_h, const handle_t& right_h) {
    uint64_t left_rank = number_bool_packing::unpack_number(left_h);
    uint64_t right_rank = number_bool_packing::unpack_number(right_h);
    auto& left_node = node_v.at(left_rank);
    auto& right_node = node_v.at(right_rank);
    bool left_rev = get_is_reverse(left_h);
    bool right_rev = get_is_reverse(right_h);
    uint64_t right_relative = edge_to_delta(right_h, left_h);
    uint64_t left_relative = edge_to_delta(left_h, right_h);

    // remove the edge from both sides
    nid_t right_node_id = get_id(right_h);
    nid_t left_node_id = get_id(left_h);

    std::vector<uint64_t> left_node_edges = left_node.edges();
    bool found_edge = false;
    for (uint64_t i = 0; i < left_node_edges.size(); ) {
        uint64_t other_id = edge_delta_to_id(left_node_id, left_node_edges.at(i++));
        uint8_t packed_edge = left_node_edges.at(i++);
        bool on_rev = edge_helper::unpack_on_rev(packed_edge);
        bool other_rev = edge_helper::unpack_other_rev(packed_edge);
        bool to_curr = edge_helper::unpack_to_curr(packed_edge);
        if (left_rev != on_rev) {
            other_rev ^= 1;
            to_curr ^= 1;
        }
        if (other_id == right_node_id && other_rev == right_rev) {
            left_node.remove_edge((i-2)/2); //convert to edge rank 
            found_edge = true;
            break;
        }
    }

    std::vector<uint64_t> right_node_edges = right_node.edges();
    for (uint64_t i = 0; i < right_node_edges.size(); ) {
        uint64_t other_id = edge_delta_to_id(right_node_id, right_node_edges.at(i++));
        uint8_t packed_edge = right_node_edges.at(i++);
        bool on_rev = edge_helper::unpack_on_rev(packed_edge);
        bool other_rev = edge_helper::unpack_other_rev(packed_edge);
        bool to_curr = edge_helper::unpack_to_curr(packed_edge);
        if (right_rev != on_rev) {
            other_rev ^= 1;
            to_curr ^= 1;
        }
        if (other_id == left_node_id && other_rev == left_rev) {
            right_node.remove_edge((i-2)/2); // convert to edge rank
            found_edge = true;
            break;
        }
    }

    _edge_count -= found_edge;
}
        
/// Remove all nodes, edges, and paths.
void ODGI::clear(void) {
    suc_bv null_bv;
    _max_node_rank = 0;
    _min_node_rank = std::numeric_limits<decltype(_min_node_rank)>::max();
    _id_increment = 1;
    _node_count = 0;
    _edge_count = 0;
    _path_count = 0;
    _path_rank_next = 0;
    deleted_node_bv = null_bv;
    // Deleted node BV needs a final delimiter.
    deleted_node_bv.push_back(1);
    node_v.clear();
    path_metadata_v.clear();
    path_name_map.clear();
}

void ODGI::clear_paths(void) {
    for_each_handle([&](const handle_t& handle) {
            node_t& node = node_v.at(number_bool_packing::unpack_number(handle));
            node.clear_path_steps();
        });
    _path_count = 0;
    _path_rank_next = 0;
    path_metadata_v.clear();
    path_name_map.clear();
}
    
/// Swap the nodes corresponding to the given handles, in the ordering used
/// by for_each_handle when looping over the graph. Other handles to the
/// nodes being swapped must not be invalidated. If a swap is made while
/// for_each_handle is running, it affects the order of the handles
/// traversed during the current traversal (so swapping an already seen
/// handle to a later handle's position will make the seen handle be visited
/// again and the later handle not be visited at all).
void ODGI::swap_handles(const handle_t& a, const handle_t& b) {
    //assert(false);
}

void ODGI::optimize(bool allow_id_reassignment) {
    apply_ordering({}, true);
}

/// Reorder the graph's internal structure to match that given.
/// Optionally compact the id space of the graph to match the ordering, from 1->|ordering|.
void ODGI::apply_ordering(const std::vector<handle_t>& order_in, bool compact_ids) {
    ODGI ordered;
    // if we're given an empty order, just compact the ids based on our ordering
    const std::vector<handle_t>* order;
    std::vector<handle_t> base_order;
    if (order_in.empty()) {
        base_order.reserve(get_node_count());
        for_each_handle([&](const handle_t& handle) {
                base_order.push_back(handle);
            });
        order = &base_order;
    } else {
        order = &order_in;
    }
    // nodes
    std::vector<nid_t> ids;
    uint64_t max_handle_rank = 0;
    uint64_t min_handle_rank = std::numeric_limits<uint64_t>::max();
    for_each_handle([&](const handle_t& handle) {
        max_handle_rank = std::max(max_handle_rank,
                                   number_bool_packing::unpack_number(handle));
        min_handle_rank = std::min(min_handle_rank,
                                   number_bool_packing::unpack_number(handle));
    });
    if (max_handle_rank > 0) {
        ids.resize(max_handle_rank - min_handle_rank + 1);
    }
    // establish id mapping
    if (compact_ids) {
        for (uint64_t i = 0; i < order->size(); ++i) {
            ids[number_bool_packing::unpack_number(order->at(i)) - min_handle_rank] = i+1;
        }
    } else {
        for (uint64_t i = 0; i < order->size(); ++i) {
            auto& handle = order->at(i);
            ids[number_bool_packing::unpack_number(handle) - min_handle_rank] = get_id(handle);
        }
    }
    // nodes
    for (auto& handle : *order) {
        ordered.create_handle(get_sequence(handle),
                              ids[number_bool_packing::unpack_number(handle) - min_handle_rank]);
    }
    // edges
    for (auto& handle : *order) {
        follow_edges(handle, false, [&](const handle_t& h) {
                ordered.create_edge(ordered.get_handle(ids[number_bool_packing::unpack_number(handle) - min_handle_rank],
                                                       get_is_reverse(handle)),
                                    ordered.get_handle(ids[number_bool_packing::unpack_number(h) - min_handle_rank],
                                                       get_is_reverse(h)));
            });
        follow_edges(handle, true, [&](const handle_t& h) {
                ordered.create_edge(ordered.get_handle(ids[number_bool_packing::unpack_number(h) - min_handle_rank],
                                                       get_is_reverse(h)),
                                    ordered.get_handle(ids[number_bool_packing::unpack_number(handle) - min_handle_rank],
                                                       get_is_reverse(handle)));
            });
    }
    // paths
    for_each_path_handle([&](const path_handle_t& old_path) {
            path_handle_t new_path = ordered.create_path_handle(get_path_name(old_path));
            for_each_step_in_path(old_path, [&](const step_handle_t& step) {
                    handle_t old_handle = get_handle_of_step(step);
                    handle_t new_handle = ordered.get_handle(ids[number_bool_packing::unpack_number(old_handle) - min_handle_rank],
                                                             get_is_reverse(old_handle));
                    ordered.append_step(new_path, new_handle);
                });
        });
    *this = ordered;
}

void ODGI::apply_path_ordering(const std::vector<path_handle_t>& order) {
    ODGI ordered;
    // copy nodes
    for_each_handle([&](const handle_t& handle) {
            ordered.create_handle(get_sequence(handle), get_id(handle));
        });
    // copy edges
    for_each_handle([&](const handle_t& handle) {
            follow_edges(handle, false, [&](const handle_t& h) {
                    ordered.create_edge(ordered.get_handle(get_id(handle), get_is_reverse(handle)),
                                        ordered.get_handle(get_id(h), get_is_reverse(h)));
                });
            follow_edges(flip(handle), false, [&](const handle_t& h) {
                    ordered.create_edge(ordered.get_handle(get_id(handle), get_is_reverse(flip(handle))),
                                        ordered.get_handle(get_id(h), get_is_reverse(h)));
                });
        });
    // add the paths in order
    for (auto& path_handle : order) {
        path_handle_t new_path = ordered.create_path_handle(get_path_name(path_handle));
        for_each_step_in_path(path_handle, [&](const step_handle_t& step) {
                handle_t old_handle = get_handle_of_step(step);
                handle_t new_handle = ordered.get_handle(get_id(old_handle), get_is_reverse(old_handle));
                ordered.append_step(new_path, new_handle);
            });
    }
    *this = ordered;
}

/// Alter the node that the given handle corresponds to so the orientation
/// indicated by the handle becomes the node's local forward orientation.
/// Rewrites all edges pointing to the node and the node's sequence to
/// reflect this. Invalidates all handles to the node (including the one
/// passed). Returns a new, valid handle to the node in its new forward
/// orientation. Note that it is possible for the node's ID to change.
/// Updates all stored paths. May change the ordering of the underlying
/// graph.
handle_t ODGI::apply_orientation(const handle_t& handle) {

    // do nothing if we're already in the right orientation
    if (!get_is_reverse(handle)) return handle;
    handle_t fwd_handle = flip(handle);
    handle_t rev_handle = handle;
    // store edges
    vector<handle_t> edges_fwd_fwd;
    vector<handle_t> edges_fwd_rev;
    vector<handle_t> edges_rev_fwd;
    vector<handle_t> edges_rev_rev;
    follow_edges(fwd_handle, false, [&](const handle_t& h) {
            edges_fwd_fwd.push_back(h);
        });
    follow_edges(fwd_handle, true, [&](const handle_t& h) {
            edges_fwd_rev.push_back(h);
        });
    follow_edges(rev_handle, false, [&](const handle_t& h) {
            edges_rev_fwd.push_back(h);
        });
    follow_edges(rev_handle, true, [&](const handle_t& h) {
            edges_rev_rev.push_back(h);
        });
    for (auto& h : edges_fwd_fwd) destroy_edge(fwd_handle, h);
    for (auto& h : edges_fwd_rev) destroy_edge(h, fwd_handle);
    for (auto& h : edges_rev_fwd) destroy_edge(rev_handle, h);
    for (auto& h : edges_rev_rev) destroy_edge(h, rev_handle);
    // save the sequence's reverse complement, which we will use to add the new handle
    //const std::string seq = get_sequence(handle);
    // we have the technology. we can rebuild it.
    // replace the handle sequence
    //set_handle_sequence(handle, seq);
    auto& node = node_v.at(number_bool_packing::unpack_number(handle));

    // flip the node sequence
    node.set_sequence(get_sequence(handle));

    // we need to flip any stored path fronts and backs as well
    std::pair<std::map<uint64_t, std::pair<uint64_t, bool>>, // path fronts
              std::map<uint64_t, std::pair<uint64_t, bool>>> // path backs
        path_rewrites = node.flip_paths(path_begin_marker, path_end_marker);
    for (auto& r : path_rewrites.first) {
        auto& p = find_metadata(rank_to_path(r.first));
        step_handle_t step;
        as_integers(step)[0]
            = as_integer(number_bool_packing::pack(number_bool_packing::unpack_number(handle),
                                                   r.second.second));
        as_integers(step)[1] = r.second.first;
        p.first = step;
    }
    for (auto& r : path_rewrites.second) {
        auto& p = find_metadata(rank_to_path(r.first));
        step_handle_t step;
        as_integers(step)[0]
            = as_integer(number_bool_packing::pack(number_bool_packing::unpack_number(handle),
                                                   r.second.second));
        as_integers(step)[1] = r.second.first;
        p.last = step;
    }
    // reconnect it to the graph
    for (auto& h : edges_fwd_fwd) {
        create_edge(handle, h);
    }
    for (auto& h : edges_fwd_rev) {
        create_edge(h, handle);
    }
    for (auto& h : edges_rev_fwd) {
        create_edge(flip(handle), h);
    }
    for (auto& h : edges_rev_rev) {
        create_edge(h, flip(handle));
    }
    return flip(handle);
}

void ODGI::set_handle_sequence(const handle_t& handle, const std::string& seq) {
    assert(seq.size());
    node_v[number_bool_packing::unpack_number(handle)].set_sequence(seq);
}
    
/// Split a handle's underlying node at the given offsets in the handle's
/// orientation. Returns all of the handles to the parts. Other handles to
/// the node being split may be invalidated. The split pieces stay in the
/// same local forward orientation as the original node, but the returned
/// handles come in the order and orientation appropriate for the handle
/// passed in.
/// Updates stored paths.
std::vector<handle_t> ODGI::divide_handle(const handle_t& handle, const std::vector<size_t>& offsets) {
    // convert the offsets to the forward strand, if needed
    std::vector<uint64_t> fwd_offsets = { 0 };
    uint64_t length = get_length(handle);
    if (get_is_reverse(handle)) {
        for (auto& o : offsets) fwd_offsets.push_back(length-o);
    } else {
        for (auto& o : offsets) fwd_offsets.push_back(o);
    }
    std::sort(fwd_offsets.begin(), fwd_offsets.end());
    handle_t fwd_handle = get_is_reverse(handle) ? flip(handle) : handle;
    // break it into the given pieces by building up the new node sequences
    std::string seq = get_sequence(fwd_handle);
    fwd_offsets.push_back(seq.size());
    std::vector<std::string> seqs;
    for (uint64_t i = 0; i < fwd_offsets.size()-1; ++i) {
        seqs.push_back(seq.substr(fwd_offsets[i], fwd_offsets[i+1]-fwd_offsets[i]));
    }
    // make the handles
    std::vector<handle_t> handles;
    for (auto& s : seqs) {
        handles.push_back(create_handle(s));
    }
    // and record their reverse, for use in path fixup
    std::vector<handle_t> rev_handles;
    for (auto& h : handles) {
        rev_handles.push_back(flip(h));
    }
    std::reverse(rev_handles.begin(), rev_handles.end());
    // connect the pieces head to tail
    for (uint64_t i = 0; i < handles.size()-1; ++i) {
        create_edge(handles[i], handles[i+1]);
    }
    // collect the handle's path context
    vector<step_handle_t> steps;
    for_each_step_on_handle(handle, [&](const step_handle_t& step) {
            steps.push_back(step);
        });
    // reverse the order to allow for safe rewriting
    std::reverse(steps.begin(), steps.end());
    // replace path steps with the new handles
    for (auto& step : steps) {
        handle_t h = get_handle_of_step(step);
        if (get_is_reverse(h)) {
            rewrite_segment(step, get_next_step(step), rev_handles);
        } else {
            rewrite_segment(step, get_next_step(step), handles);
        }
    }
    node_v.at(number_bool_packing::unpack_number(handle)).clear_path_steps();
    // collect the context of the handle
    vector<handle_t> edges_fwd_fwd;
    vector<handle_t> edges_fwd_rev;
    vector<handle_t> edges_rev_fwd;
    vector<handle_t> edges_rev_rev;
    follow_edges(fwd_handle, false, [&](const handle_t& h) {
            edges_fwd_fwd.push_back(h);
        });
    follow_edges(fwd_handle, true, [&](const handle_t& h) {
            edges_fwd_rev.push_back(h);
        });
    follow_edges(flip(fwd_handle), false, [&](const handle_t& h) {
            edges_rev_fwd.push_back(h);
        });
    follow_edges(flip(fwd_handle), true, [&](const handle_t& h) {
            edges_rev_rev.push_back(h);
        });
    // destroy the handle
    destroy_handle(fwd_handle);
    // connect the ends to the previous context
    for (auto& h : edges_fwd_fwd) create_edge(handles.back(), h);
    for (auto& h : edges_fwd_rev) create_edge(h, handles.front());
    for (auto& h : edges_rev_fwd) create_edge(rev_handles.back(), h);
    for (auto& h : edges_rev_rev) create_edge(h, rev_handles.front());
    return get_is_reverse(handle) ? rev_handles : handles;
}

handle_t ODGI::combine_handles(const std::vector<handle_t>& handles) {
    std::string seq;
    for (auto& handle : handles) {
        seq.append(get_sequence(handle));
    }
    handle_t combined = create_handle(seq);
    // relink the inbound and outbound nodes
    // get the edge context
    vector<handle_t> edges_fwd_fwd;
    vector<handle_t> edges_fwd_rev;
    vector<handle_t> edges_rev_fwd;
    vector<handle_t> edges_rev_rev;
    follow_edges(handles.back(), false, [&](const handle_t& h) {
            edges_fwd_fwd.push_back(h);
        });
    follow_edges(handles.front(), true, [&](const handle_t& h) {
            edges_fwd_rev.push_back(h);
        });
    // destroy the old handles
    for (auto& handle : handles) {
        destroy_handle(handle);
    }
    // connect the ends to the previous context
    // check that we're not trying to make edges that connect back with the nodes in the component
    // there are three cases
    // self looping, front and rear inverting
    for (auto& h : edges_fwd_fwd) {
        if (h == handles.front()) {
            create_edge(combined, combined);
        } else if (h == flip(handles.back())) {
            create_edge(combined, flip(combined));
        } else {
            create_edge(combined, h);
        }
    }
    for (auto& h : edges_fwd_rev) {
        if (h == handles.back()) {
            create_edge(combined, combined);
        } else if (h == flip(handles.front())) {
            create_edge(flip(combined), combined);
        } else {
            create_edge(h, combined);
        }
    }
    return combined;
}

/*
 * This is the interface for a handle graph with embedded paths where the paths can be modified.
 * Note that if the *graph* can also be modified, the implementation will also
 * need to inherit from MutableHandleGraph, via the combination
 * MutablePathMutableHandleGraph interface.
 * TODO: This is a very limited interface at the moment. It will probably need to be extended.
 */
    
/**
 * Destroy the given path. Invalidates handles to the path and its node steps.
 */
void ODGI::destroy_path(const path_handle_t& path) {
    // select everything with that handle in the path_handle_wt
    std::vector<step_handle_t> path_v;
    for_each_step_in_path(path, [this,&path_v](const step_handle_t& step) {
            path_v.push_back(step);
        });
    // this is order dependent...
    // we need to destroy steps in their reverse ranks
    std::sort(path_v.begin(), path_v.end(), [](const step_handle_t& a, const step_handle_t& b) {
            return as_integers(a)[0] == as_integers(b)[0] && as_integers(a)[1] > as_integers(b)[1]
                || as_integers(a)[0] < as_integers(b)[0]; });
    for (auto& step : path_v) {
        destroy_step(step);
    }
    // erase it from the name index if destroy_step didn't get it.
    auto f = path_name_map.find(get_path_name(path));
    if (f != path_name_map.end()) {
        path_name_map.erase(f);
    }
    auto& p = find_metadata(path);
    // our length should be 0
    assert(p.length == 0);
    --_path_count;
}

/**
 * Create a path with the given name. The caller must ensure that no path
 * with the given name exists already, or the behavior is undefined.
 * Returns a handle to the created empty path. Handles to other paths must
 * remain valid.
 */
path_handle_t ODGI::create_path_handle(const std::string& name, bool is_circular) {
    path_handle_t path = rank_to_path(_path_rank_next++);
    path_name_map[name] = path;
    path_metadata_v.emplace_back();
    // Make sure we are putting the metadata where we will find it.
    assert(path_metadata_v.size() == _path_rank_next);
    auto& p = find_metadata(path);
    // Path metadata starts out with 0 length and zeroed step handles. 
    p.name = name;
    p.is_circular = is_circular;
    ++_path_count;
    return path;
}

step_handle_t ODGI::create_step(const path_handle_t& path, const handle_t& handle) {
    // where are we going to insert?
    uint64_t rank_on_handle = get_step_count(handle);
    // build our step
    step_handle_t step;
    as_integers(step)[0] = as_integer(handle);
    as_integers(step)[1] = rank_on_handle;
    auto& node = node_v.at(number_bool_packing::unpack_number(handle));
    node.add_path_step(path_to_rank(path), get_is_reverse(handle),
                       path_begin_marker, 0, path_end_marker, 0);
    return step;
}

void ODGI::link_steps(const step_handle_t& from, const step_handle_t& to) {
    path_handle_t path = get_path_handle_of_step(from);
    assert(path == get_path_handle_of_step(to));
    const handle_t& from_handle = get_handle_of_step(from);
    const handle_t& to_handle = get_handle_of_step(to);
    const uint64_t& from_rank = as_integers(from)[1];
    const uint64_t& to_rank = as_integers(to)[1];
    node_t& from_node = node_v.at(number_bool_packing::unpack_number(from_handle));
    node_t::step_t from_step = from_node.get_path_step(from_rank);
    from_step.set_next_id(edge_to_delta(from_handle, to_handle) + (path_end_marker + 1));
    from_step.set_next_rank(to_rank);
    from_node.set_path_step(from_rank, from_step);
    node_t& to_node = node_v.at(number_bool_packing::unpack_number(to_handle));
    node_t::step_t to_step = to_node.get_path_step(to_rank);
    to_step.set_prev_id(edge_to_delta(to_handle, from_handle) + (path_end_marker + 1));
    to_step.set_prev_rank(from_rank);
    to_node.set_path_step(to_rank, to_step);
}

void ODGI::destroy_step(const step_handle_t& step_handle, bool clean_up_empty_path) {
    // erase reference to this step
    bool has_prev = has_linear_previous_step(step_handle);
    bool has_next = has_linear_next_step(step_handle);
    if (!has_prev && !has_next && clean_up_empty_path) {
        // we're about to erase the path, so we need to clean up the path metadata record
        path_handle_t path = get_path_handle_of_step(step_handle);
        path_name_map.erase(get_path_name(path));
        find_metadata(path) = path_metadata_t();
    } else {
        if (has_prev) {
            auto step = get_previous_step(step_handle);
            node_t& step_node = node_v.at(number_bool_packing::unpack_number(get_handle_of_step(step)));
            uint64_t step_rank = as_integers(step)[1];
            node_t::step_t node_step = step_node.get_path_step(step_rank);
            //std::cerr << "destroy prev links " << step_node.id() << std::endl;
            node_step.set_next_id(path_end_marker);
            node_step.set_next_rank(0);
            step_node.set_path_step(step_rank, node_step);
        } else if (has_next) {
            auto step = get_next_step(step_handle);
            auto& p = find_metadata(get_path_handle_of_step(step));
            p.first = step;
        }
        if (has_next) {
            auto step = get_next_step(step_handle);
            node_t& step_node = node_v.at(number_bool_packing::unpack_number(get_handle_of_step(step)));
            uint64_t step_rank = as_integers(step)[1];
            node_t::step_t node_step = step_node.get_path_step(step_rank);
            //std::cerr << "destroy next links " << step_node.id() << std::endl;
            node_step.set_prev_id(path_begin_marker);
            node_step.set_prev_rank(0);
            step_node.set_path_step(step_rank, node_step);
        } else if (has_prev) {
            auto step = get_previous_step(step_handle);
            auto& p = find_metadata(get_path_handle_of_step(step));
            p.last = step;
        }
    }
    // update other records on this path on this node
    handle_t handle = get_handle_of_step(step_handle);
    bool seen_curr = false;
    for_each_step_on_handle(handle, [&](const step_handle_t& step) {
            if (seen_curr) {
                decrement_rank(step);
            }
            if (step == step_handle) {
                seen_curr = true;
            }
        });
    node_t& curr_node = node_v.at(number_bool_packing::unpack_number(get_handle_of_step(step_handle)));
    curr_node.remove_path_step(as_integers(step_handle)[1]);
    
    if (curr_node.path_count() == 0 && graph_id_hidden_set.count(get_id(handle))) {
        // We have freed up a hidden handle from its last reference. Clean it
        // up and get rid of the now-useless node.
        destroy_handle(handle);
    }
}

step_handle_t ODGI::prepend_step(const path_handle_t& path, const handle_t& to_append) {
    // get the last step
    auto& p = find_metadata(path);
    // create the new step
    step_handle_t new_step = create_step(path, to_append);
    if (!p.length) {
        p.last = new_step;
    } else {
        step_handle_t first_step = path_begin(path);
        // link it to the last step
        link_steps(new_step, first_step);
    }
    // point to the new last step
    p.first = new_step;
    // update our step count
    ++p.length;
    return new_step;
}

step_handle_t ODGI::append_step(const path_handle_t& path, const handle_t& to_append) {
    // get the last step
    auto& p = find_metadata(path);
    // create the new step
    step_handle_t new_step = create_step(path, to_append);
    if (!p.length) {
        p.first = new_step;
    } else {
        step_handle_t last_step = path_back(path);
        // link it to the last step
        link_steps(last_step, new_step);
    }
    // point to the new last step
    p.last = new_step;
    // update our step count
    ++p.length;
    return new_step;
}

/// helper to handle the case where we remove an step from a given path
/// on a node that has other steps from the same path, thus invalidating the
/// ranks used to refer to it
void ODGI::decrement_rank(const step_handle_t& step_handle) {
    // what is the actual rank of this step?
    //std::cerr << "in decrement rank " << get_handle_of_step(step_handle) << ":" << as_integers(step_handle)[1] << std::endl;
    if (has_linear_previous_step(step_handle)) {
        auto step = get_previous_step(step_handle);
        // decrement the rank information
        node_t& step_node = node_v.at(number_bool_packing::unpack_number(get_handle_of_step(step)));
        uint64_t step_rank = as_integers(step)[1];
        node_t::step_t node_step = step_node.get_path_step(step_rank);
        node_step.set_next_rank(node_step.next_rank()-1);
        step_node.set_path_step(step_rank, node_step);
    } else {
        // update path metadata
        auto& p = find_metadata(get_path_handle_of_step(step_handle));
        --as_integers(p.first)[1];
    }
    if (has_linear_next_step(step_handle)) {
        auto step = get_next_step(step_handle);
        node_t& step_node = node_v.at(number_bool_packing::unpack_number(get_handle_of_step(step)));
        uint64_t step_rank = as_integers(step)[1];
        node_t::step_t node_step = step_node.get_path_step(step_rank);
        node_step.set_prev_rank(node_step.prev_rank()-1);
        step_node.set_path_step(step_rank, node_step);
    } else {
        // update path metadata
        auto& p = find_metadata(get_path_handle_of_step(step_handle));
        --as_integers(p.last)[1];
    }
}

/// Modify the given step handle to point to the given handle.
void ODGI::set_handle_of_step(step_handle_t& step_handle, const handle_t& handle) const {
    as_integers(step_handle)[0] = as_integer(handle);
}

/// Find the path metadata for a path
ODGI::path_metadata_t& ODGI::find_metadata(const path_handle_t& path) {
    // Nobody but us should touch path_metadata_v indexing!
    return path_metadata_v.at(path_to_rank(path));
}

/// Find the metadata for a path, read-only
const ODGI::path_metadata_t& ODGI::find_metadata(const path_handle_t& path) const {
    // Nobody but us should touch path_metadata_v indexing!
    return path_metadata_v.at(path_to_rank(path));
}

/// Get the 0-based rank encoded by a path handle
size_t ODGI::path_to_rank(const path_handle_t& path) const {
    // Nobody but us should unpack path handles!
    return as_integer(path) - 1;
}

/// Get the path handle encodign the given 0-based rank
path_handle_t ODGI::rank_to_path(size_t rank) const {
    // Nobody but us should pack path handles!
    return as_path_handle(rank + 1);
}

/// Add the offset to the number packed in the given handle, and return a
/// new modified handle.
handle_t ODGI::add_to_number(const handle_t& handle, int64_t offset) const {
    return number_bool_packing::pack(number_bool_packing::unpack_number(handle) + offset,
        number_bool_packing::unpack_bit(handle));
}

// Insert a visit to a node to the given path between the given steps.
step_handle_t ODGI::insert_step(const step_handle_t& before, const step_handle_t& after, const handle_t& to_insert) {
    auto p = rewrite_segment(after, after, { to_insert });
    return get_next_step(p.first);
}

/// reassign the given step to the new handle
step_handle_t ODGI::set_step(const step_handle_t& step_handle, const handle_t& assign_to) {
    return rewrite_segment(step_handle, get_next_step(step_handle), { assign_to }).first;
}

/// Replace the path range with the new segment
std::pair<step_handle_t, step_handle_t> ODGI::rewrite_segment(const step_handle_t& segment_begin,
                                                                 const step_handle_t& segment_end,
                                                                 const std::vector<handle_t>& new_segment) {
    
    // Detect if we are crossing the weld on a circular path
    bool crossed_endpoint = false;
    
    // collect the steps to replace
    std::vector<step_handle_t> steps;
    for (step_handle_t step = segment_begin; step != segment_end; step = get_next_step(step)) {
        steps.push_back(step);
        if (step != segment_begin && !has_linear_previous_step(step)) {
            // We crossed over the end of the path. Path ought to be circular,
            // and we need to take care when hooking things up not to make a
            // circular data structure.
            crossed_endpoint = true;
        }
    }
   
    // get the path metadata
    path_handle_t path = get_path_handle_of_step(segment_begin);
    auto& path_meta = find_metadata(path);
    
    // find the before and after steps, which we'll link into, if they exist
    // Ignore circularity; that is handled outside our link-based data structure.
    bool has_predecessor = (path_meta.length != 0 && has_linear_previous_step(segment_begin));
    bool has_successor = !is_path_end(segment_end);
    step_handle_t before;
    step_handle_t after;
    if (has_predecessor) {
        before = get_previous_step(segment_begin);
    }
    if (has_successor) {
        after = segment_end;
    }
    
    // sort the steps so that we destroy from higher to lower ranks
    std::sort(steps.begin(), steps.end(), [](const step_handle_t& a, const step_handle_t& b) {
            return as_integers(a)[0] == as_integers(b)[0] && as_integers(a)[1] > as_integers(b)[1]
                || as_integers(a)[0] < as_integers(b)[0]; });
      
    // delete the previous steps, but don't clean up the whole path if it becomes empty
    for (auto& step : steps) {
        destroy_step(step, false);
    }
    path_meta.length -= steps.size();
    // create the new steps
    std::vector<step_handle_t> new_steps;
    for (auto& handle : new_segment) {
        new_steps.push_back(create_step(path, handle));
    }
    path_meta.length += new_steps.size();
    
    if (new_steps.size()) {
        // link new steps together
        for (uint64_t i = 0; i + 1 < new_steps.size(); ++i) {
            link_steps(new_steps[i], new_steps[i+1]);
        }
        if (has_predecessor && !crossed_endpoint) {
            // A step exists before what we rewrote away, and we aren't
            // rewriting across a circularization point and becoming the new
            // beginning of the path.
            link_steps(before, new_steps.front());
        } else {
            path_meta.first = new_steps.front();
        }
        if (has_successor) {
            link_steps(new_steps.back(), after);
        } else {
            path_meta.last = new_steps.back();
        }
        // Make sure we wired up to the end, or the end is the path's end and we get it correctly.
        assert(get_next_step(new_steps.back()) == segment_end);
        return make_pair(new_steps.front(), segment_end);
    } else {
        // An empty range is two of the past-end step we started with.
        return make_pair(segment_end, segment_end);
    }
}

void ODGI::display(void) const {
    std::cerr << "------ graph state ------" << std::endl;

    std::cerr << "_id_increment = " << _id_increment << std::endl;
    std::cerr << "_max_node_rank = " << _max_node_rank << std::endl;
    std::cerr << "_min_node_rank = " << _min_node_rank << std::endl;

    //std::cerr << "graph_id_map" << "\t";
    //for (auto& k : graph_id_map) std::cerr << k.first << "->" << k.second << " "; std::cerr << std::endl;
    std::cerr << "node_v" << "\t";
    for (uint64_t i = 0; i < node_v.size(); ++i) {
        auto& node = node_v.at(i);
        nid_t node_id = rank_to_id(i);
        std::cerr << node_id << ":" << node.sequence() << " ";
        const std::vector<uint64_t> node_edges = node.edges();
        for (uint64_t j = 0; j < node_edges.size(); ++j) {
            std::cerr << node_edges.at(j) << ",";
        }
        std::cerr << " _ ";
        const std::vector<node_t::step_t> steps = node.get_path_steps();
        for (auto& step : steps) {
            std::cerr << step.path_id() << ":"
                      << step.is_rev() << ":"
                      << (step.prev_id() == path_begin_marker ? "#" : std::to_string(edge_delta_to_id(node_id, step.prev_id() - (path_end_marker + 1)))) << ":"
                      << step.prev_rank() << ":"
                      << (step.next_id() == path_end_marker ? "$" : std::to_string(edge_delta_to_id(node_id, step.next_id() - (path_end_marker + 1)))) << ":"
                      << step.next_rank() << " ";
        }
        std::cerr << " | ";
    }
    std::cerr << std::endl;
    std::cerr << "deleted_node_bv" << "\t";
    for (uint64_t i = 0; i < deleted_node_bv.size(); ++i) std::cerr << deleted_node_bv.at(i) << " "; std::cerr << std::endl;
    /// Ordered across the nodes in graph_id_iv, stores the path ids (1-based) at each
    /// segment in seq_wt, delimited by 0, one for each path step (node traversal).
    std::cerr << "path_metadata" << "\t";
    for (uint64_t q = 0; q < path_metadata_v.size(); ++q) {
        auto& p = path_metadata_v.at(q);
        std::cerr << q << ":" << p.name << ":"
                  << as_integers(p.first)[0] << "(" <<  get_id(as_handle(as_integers(p.first)[0])) << ")/" << as_integers(p.first)[1] << "->"
                  << as_integers(p.last)[0] << "(" <<  get_id(as_handle(as_integers(p.last)[0])) << ")/" << as_integers(p.last)[1] << " ";
    } std::cerr << std::endl;

    // not dumped...
    /// Stores path names in their internal order, delimited by '$'
    //dyn::wt_fmi path_name_fmi;
    /// Marks the beginning of each path name
    //dyn::suc_bv path_name_bv;

}

void ODGI::to_gfa(std::ostream& out) const {
    out << "H\tVN:Z:1.0" << std::endl;
    // for each node
    for_each_handle([&out,this](const handle_t& h) {
            out << "S\t" << get_id(h) << "\t" << get_sequence(h) << std::endl;
            {
                // use this direct iteration to avoid double counting edges
                // we only consider write the edges relative to their start
                const node_t& node = node_v.at(number_bool_packing::unpack_number(h));
                bool is_rev = get_is_reverse(h);
                nid_t node_id = get_id(h);
                const std::vector<uint64_t> node_edges = node.edges();
                for (uint64_t i = 0; i < node_edges.size(); i+=2) {
                    // unpack the edge
                    uint64_t other_id = edge_delta_to_id(node_id, node_edges.at(i));
                    uint8_t packed_edge = node_edges.at(i+1);
                    bool on_rev = edge_helper::unpack_on_rev(packed_edge);
                    bool other_rev = edge_helper::unpack_other_rev(packed_edge);
                    bool to_curr = edge_helper::unpack_to_curr(packed_edge);
                    if (!to_curr) {
                        out << "L\t" << node_id << "\t"
                            << (on_rev?"-":"+")
                            << "\t" << other_id << "\t"
                            << (other_rev?"-":"+")
                            << "\t0M" << std::endl;
                    }
                }
            }
        });
    for_each_path_handle([&out,this](const path_handle_t& p) {
            //step_handle_t step = path_begin(p);
            out << "P\t" << get_path_name(p) << "\t";
            for_each_step_in_path(p, [this,&out](const step_handle_t& step) {
                    handle_t h = get_handle_of_step(step);
                    out << get_id(h) << (get_is_reverse(h)?"-":"+");
                    if (has_next_step(step)) out << ",";
                });
            out << "\t";
            for_each_step_in_path(p, [this,&out](const step_handle_t& step) {
                    out << get_length(get_handle_of_step(step)) << "M";
                    if (has_next_step(step)) out << ",";
                });
            out << std::endl;
        });
}

void ODGI::to_gfa(const string& filename) const {
    if (filename == "-") {
        to_gfa(std::cout);
    } else {
        ofstream out(filename);
        if (!out) {
            throw runtime_error("Could not open " + filename + " for writing");
        }
        to_gfa(out);
        if (!out) {
            throw runtime_error("Could not write GFA to " + filename);
        }
    }
}

string ODGI::to_gfa() const {
    stringstream ss;
    to_gfa(ss);
    return ss.str();
}

long long int ODGI::serialize_and_measure(std::ostream& out) const {
    //rebuild_id_handle_mapping();
    // TODO: every write here needs an endianness-converting function!
    long long int written = 0;
    out.write((char*)&_max_node_rank,sizeof(_max_node_rank));
    written += sizeof(_max_node_rank);
    out.write((char*)&_min_node_rank,sizeof(_min_node_rank));
    written += sizeof(_min_node_rank);
    out.write((char*)&_node_count,sizeof(_node_count));
    written += sizeof(_node_count);
    out.write((char*)&_edge_count,sizeof(_edge_count));
    written += sizeof(_edge_count);
    out.write((char*)&_path_count,sizeof(_path_count));
    written += sizeof(_path_count);
    out.write((char*)&_path_rank_next,sizeof(_path_rank_next));
    written += sizeof(_path_rank_next);
    out.write((char*)&_deleted_node_count,sizeof(_deleted_node_count));
    written += sizeof(_deleted_node_count);
    out.write((char*)&_id_increment,sizeof(_id_increment));
    written += sizeof(_id_increment);
    size_t hidden_count = graph_id_hidden_set.size();
    out.write((char*)&hidden_count,sizeof(hidden_count));
    written += sizeof(hidden_count);
    for (const nid_t& hidden : graph_id_hidden_set) {
        out.write((const char*)&hidden,sizeof(hidden));
        written += sizeof(hidden);
    }
    assert(_node_count == node_v.size());
    for (auto& node : node_v) {
        written += node.serialize(out);
    }
    written += deleted_node_bv.serialize(out);
    size_t i = path_metadata_v.size();
    out.write((char*)&i,sizeof(size_t));
    written += sizeof(size_t);
    for (auto m : path_metadata_v) {
        out.write((char*)&m.length,sizeof(m.length));
        written += sizeof(m.length);
        out.write((char*)&m.first,sizeof(m.first));
        written += sizeof(m.first);
        out.write((char*)&m.last,sizeof(m.last));
        written += sizeof(m.last);
        i = m.name.size();
        out.write((char*)&i,sizeof(size_t));
        written += sizeof(size_t);
        out.write((char*)m.name.c_str(),m.name.size());
        written += m.name.size();
    }
    i = path_name_map.size();
    out.write((char*)&i,sizeof(size_t));
    written += sizeof(size_t);
    for (auto& p : path_name_map) {
        i = p.first.size();
        out.write((char*)&i,sizeof(size_t));
        written += sizeof(size_t);
        out.write(p.first.c_str(),p.first.size());
        written += p.first.size();
        out.write((char*)&p.second,sizeof(p.second));
        written += sizeof(p.second);
    }
    return written;
}

void ODGI::load(std::istream& in) {
    // TODO: every read here needs an endianness-converting function!
    in.read((char*)&_max_node_rank,sizeof(_max_node_rank));
    in.read((char*)&_min_node_rank,sizeof(_min_node_rank));
    in.read((char*)&_node_count,sizeof(_node_count));
    in.read((char*)&_edge_count,sizeof(_edge_count));
    in.read((char*)&_path_count,sizeof(_path_count));
    in.read((char*)&_path_rank_next,sizeof(_path_rank_next));
    in.read((char*)&_deleted_node_count,sizeof(_deleted_node_count));
    in.read((char*)&_id_increment,sizeof(_id_increment));
    size_t hidden_count;
    in.read((char*)&hidden_count,sizeof(hidden_count));
    for (size_t i = 0; i < hidden_count; i++) {
        nid_t hidden;
        in.read((char*)&hidden,sizeof(hidden));
        graph_id_hidden_set.insert(hidden);
    }
    node_v.resize(_node_count);
    for (size_t i = 0; i < _node_count; ++i) {
        auto& node = node_v[i];
        node.load(in);
    }
    deleted_node_bv.load(in);
    size_t i = 0;
    in.read((char*)&i,sizeof(size_t));
    path_metadata_v.reserve(i);
    for (size_t j = 0; j < i; ++j) {
        path_metadata_v.emplace_back();
        auto& m = path_metadata_v.back();
        in.read((char*)&m.length,sizeof(m.length));
        in.read((char*)&m.first,sizeof(m.first));
        in.read((char*)&m.last,sizeof(m.last));
        uint64_t s;
        in.read((char*)&s,sizeof(size_t));
        char n[s+1]; n[s] = '\0';
        in.read(n,s);
        m.name = string(n);
    }
    in.read((char*)&i,sizeof(size_t));
    path_name_map.reserve(i);
    for (size_t j = 0; j < i; ++j) {
        uint64_t s;
        in.read((char*)&s,sizeof(size_t));
        char k[s+1]; k[s] = '\0';
        in.read(k,s);
        path_handle_t v;
        in.read((char*)&v,sizeof(v));
        path_name_map[string(k)] = v;
    }
    
    if (_id_increment == 0) {
        // Make sure _id_increment isn't 0. Otherwise, the rank-0 unused slot will
        // try to have ID 0 and we can't use it. This upsets the create_handle
        // logic that tries to re-use the first empty slot.
        set_id_increment(1);
    }
    
}

}
