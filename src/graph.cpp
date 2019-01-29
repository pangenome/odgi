//
//  graph.cpp
//  

#include "graph.hpp"

namespace dankgraph {

/// Look up the handle for the node with the given ID in the given orientation
handle_t graph_t::get_handle(const id_t& node_id, bool is_reverse) const {
    //return handle_helper::pack(graph_id_wt.select(0, node_id), is_reverse);
    auto f = graph_id_map.find(node_id);
    assert(f != graph_id_map.end());
    return handle_helper::pack(f->second, is_reverse);
}
    
/// Get the ID from a handle
id_t graph_t::get_id(const handle_t& handle) const {
    return graph_id_pv.at(handle_helper::unpack_number(handle));
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
    return seq_bv.select1(offset+1) - seq_bv.select1(offset);
}
    
/// Get the sequence of a node, presented in the handle's local forward orientation.
std::string graph_t::get_sequence(const handle_t& handle) const {
    std::string seq;
    uint64_t offset = handle_helper::unpack_number(handle);
    for (uint64_t i = seq_bv.select1(offset); ; ++i) {
        if (seq.size() && seq_bv.at(i)) break;
        seq += int_as_dna(seq_pv.at(i));
    }
    return seq;
}
    
/// Loop over all the handles to next/previous (right/left) nodes. Passes
/// them to a callback which returns false to stop iterating and true to
/// continue. Returns true if we finished and false if we stopped early.
bool graph_t::follow_edges(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const {
    bool result = true;
    uint64_t offset = handle_helper::unpack_number(handle);
    bool is_rev = handle_helper::unpack_bit(handle);
    // NB edges are stored in canonical orientation, forward to reverse prefered
    if (!go_left && !is_rev || go_left && is_rev) {
        uint64_t edges_begin = edge_fwd_bv.select1(offset)+1;
        for (uint64_t i = edges_begin; ; ++i) {
            uint64_t x = edge_fwd_iv.at(i);
            if (x==0) break; // end of record
            uint64_t id = edge_delta_to_id(get_id(handle), x);
            bool inv = edge_fwd_inv_bv.at(i);
            handle_t handle = get_handle(id, (inv ? !is_rev : is_rev));
            result &= iteratee(handle);
            if (!result) break;
        }
    } else {
        assert(go_left && !is_rev || !go_left && is_rev);
        uint64_t edges_begin = edge_rev_bv.select1(offset)+1;
        for (uint64_t i = edges_begin; ; ++i) {
            uint64_t x = edge_rev_iv.at(i);
            if (x==0) break; // end of record
            uint64_t id = edge_delta_to_id(get_id(handle), x);
            bool inv = edge_rev_inv_bv.at(i);
            handle_t handle = get_handle(id, (inv ? !is_rev : is_rev));
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
        for (uint64_t i = 0; i < graph_id_pv.size(); ++i) {
            if (flag) continue;
            bool result = iteratee(handle_helper::pack(i,false));
#pragma omp atomic
            flag &= result;
        }
    } else {
        for (uint64_t i = 0; i < graph_id_pv.size(); ++i) {
            if (!iteratee(handle_helper::pack(i,false))) break;
        }
    }
}

void graph_t::for_each_edge(const std::function<bool(const edge_t&)>& iteratee, bool parallel){
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
    
/// Return the number of nodes in the graph
/// TODO: can't be node_count because XG has a field named node_count.
size_t graph_t::node_size(void) const {
    return graph_id_pv.size();
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
        return edge_fwd_bv.select1(offset+1) - edge_fwd_bv.select1(offset);
    } else {
        return edge_rev_bv.select1(offset+1) - edge_rev_bv.select1(offset);
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

/// A pair of handles can be used as an edge. When so used, the handles have a
/// canonical order and orientation.
edge_t graph_t::edge_handle(const handle_t& left, const handle_t& right) const {
    return make_pair(left, right);
}
    
/// Such a pair can be viewed from either inward end handle and produce the
/// outward handle you would arrive at.
/*
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
    auto f = path_name_map.find(path_name);
    if (f == path_name_map.end()) return false;
    else return true;
}
    
/// Look up the path handle for the given path name.
/// The path with that name must exist.
path_handle_t graph_t::get_path_handle(const std::string& path_name) const {
    auto f = path_name_map.find(path_name);
    assert(f != path_name_map.end());
    return as_path_handle(f->second);
}

/// Look up the name of a path from a handle to it
std::string graph_t::get_path_name(const path_handle_t& path_handle) const {
    auto f = path_metadata_map.find(as_integer(path_handle));
    assert(f != path_metadata_map.end());
    return f->second.name;
}
    
/// Returns the number of node occurrences in the path
size_t graph_t::get_occurrence_count(const path_handle_t& path_handle) const {
    auto f = path_metadata_map.find(as_integer(path_handle));
    if (f == path_metadata_map.end()) return 0;
    else return f->second.length;
}

/// Returns the number of paths stored in the graph
size_t graph_t::get_path_count(void) const {
    return _path_count;
}
    
/// Execute a function on each path in the graph
// TODO: allow stopping early?
void graph_t::for_each_path_handle(const std::function<void(const path_handle_t&)>& iteratee) const {
    for (uint64_t i = 0; i < _path_handle_next; ++i) {
        path_handle_t path = as_path_handle(i);
        if (get_occurrence_count(path) > 0) {
            iteratee(path);
        }
    }
}

void graph_t::for_each_occurrence_on_handle(const handle_t& handle, const std::function<void(const occurrence_handle_t&)>& iteratee) const {
    uint64_t handle_rank = handle_helper::unpack_number(handle);
    uint64_t begin = path_handle_wt.select(handle_rank, 0);
    uint64_t end = path_handle_wt.select(handle_rank+1, 0);
    hash_map<uint64_t, uint64_t> path_count;
    for (uint64_t i = 0; i < end-begin; ++i) {
        occurrence_handle_t occ;
        as_integers(occ)[0] = handle_rank;
        as_integers(occ)[1] = i;
        iteratee(occ);
    }
}

size_t graph_t::get_occurrence_count(const handle_t& handle) const {
    uint64_t handle_rank = handle_helper::unpack_number(handle);
    uint64_t begin = path_handle_wt.select(handle_rank, 0)+1;
    uint64_t end = path_handle_wt.select(handle_rank+1, 0);
    return end - begin;
}

uint64_t graph_t::occurrence_rank(const occurrence_handle_t& occurrence_handle) const {
    uint64_t i = as_integers(occurrence_handle)[0];
    uint64_t j = as_integers(occurrence_handle)[1];
    return path_handle_wt.select(i, 0)+1 + j;
}

/// Get a node handle (node ID and orientation) from a handle to an occurrence on a path
handle_t graph_t::get_occurrence(const occurrence_handle_t& occurrence_handle) const {
    uint64_t i = occurrence_rank(occurrence_handle);
    return handle_helper::pack(as_integers(occurrence_handle)[0], path_rev_iv.at(i));
}

/// Get a path handle (path ID) from a handle to an occurrence on a path
path_handle_t graph_t::get_path(const occurrence_handle_t& occurrence_handle) const {
    return as_path_handle(path_handle_wt.at(occurrence_rank(occurrence_handle))-1);
}

/// Get a handle to the first occurrence in a path.
/// The path MUST be nonempty.
occurrence_handle_t graph_t::get_first_occurrence(const path_handle_t& path_handle) const {
    return path_metadata_map.at(as_integer(path_handle)).first;
}
    
/// Get a handle to the last occurrence in a path
/// The path MUST be nonempty.
occurrence_handle_t graph_t::get_last_occurrence(const path_handle_t& path_handle) const {
    return path_metadata_map.at(as_integer(path_handle)).last;
}
    
/// Returns true if the occurrence is not the last occurence on the path, else false
bool graph_t::has_next_occurrence(const occurrence_handle_t& occurrence_handle) const {
    return path_next_id_iv.at(occurrence_rank(occurrence_handle)) != path_end_marker;
}
    
/// Returns true if the occurrence is not the first occurence on the path, else false
bool graph_t::has_previous_occurrence(const occurrence_handle_t& occurrence_handle) const {
    return path_next_id_iv.at(occurrence_rank(occurrence_handle)) != path_begin_marker;
}

/// Returns a handle to the next occurrence on the path, which must exist
occurrence_handle_t graph_t::get_next_occurrence(const occurrence_handle_t& occurrence_handle) const {
    uint64_t i = occurrence_rank(occurrence_handle);
    id_t curr_id = get_id(handle_helper::pack(as_integers(occurrence_handle)[0], false));
    occurrence_handle_t occ;
    //as_integers(occ)[0] = edge_delta_to_id(curr_id, path_next_id_iv.at(i));
    as_integers(occ)[0] = handle_helper::unpack_number(get_handle(edge_delta_to_id(curr_id, path_next_id_iv.at(i)-2), false));
    as_integers(occ)[1] = path_next_rank_iv.at(i);
    return occ;
}

/// Returns a handle to the previous occurrence on the path
occurrence_handle_t graph_t::get_previous_occurrence(const occurrence_handle_t& occurrence_handle) const {
    uint64_t i = occurrence_rank(occurrence_handle);
    id_t curr_id = get_id(handle_helper::pack(as_integers(occurrence_handle)[0], false));
    occurrence_handle_t occ;
    //as_integers(occ)[0] = get_handle(edge_delta_to_id(curr_id, path_prev_id_iv.at(i)), false);
    as_integers(occ)[0] = handle_helper::unpack_number(get_handle(edge_delta_to_id(curr_id, path_prev_id_iv.at(i)-2), false));
    as_integers(occ)[1] = path_prev_rank_iv.at(i);
    return occ;
}
    
////////////////////////////////////////////////////////////////////////////
// Additional optional interface with a default implementation
////////////////////////////////////////////////////////////////////////////

/// Returns true if the given path is empty, and false otherwise
bool graph_t::is_empty(const path_handle_t& path_handle) const {
    return get_occurrence_count(path_handle) == 0;
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
    // get node id as max+1
    return create_handle(sequence, _max_node_id+1);
}

handle_t graph_t::create_hidden_handle(const std::string& sequence) {
    // get node id as max+1
    uint64_t id = _max_node_id+1;
    graph_id_hidden_set.insert(id);
    ++_hidden_count;
    return create_handle(sequence, id);
}

/// Create a new node with the given id and sequence, then return the handle.
handle_t graph_t::create_handle(const std::string& sequence, const id_t& id) {
    assert(!graph_id_map.find(id) != graph_id_map.end());
    assert(id > 0);
    id_t new_id = id;
    // set new max
    _max_node_id = max(new_id, _max_node_id);
    _min_node_id = max((uint64_t)1, (uint64_t)min(new_id, _min_node_id));
    graph_id_map[new_id] = graph_id_pv.size();
    // add to graph_id_pv
    graph_id_pv.push_back(new_id);
    // append to seq_wt, delimit by 0
    for (auto c : sequence) {
        seq_pv.push_back(dna_as_int(c));
    }
    // update seq_bv
    for (uint64_t i = 0; i < sequence.size()-1; ++i) {
            seq_bv.push_back(0);
    }
    seq_bv.push_back(1); // end delimiter
    // set up delemiters for edges, for later filling
    edge_fwd_iv.push_back(0);
    edge_fwd_bv.push_back(1);
    edge_fwd_inv_bv.push_back(0);
    edge_rev_iv.push_back(0);
    edge_rev_bv.push_back(1);
    edge_rev_inv_bv.push_back(0);
    // set up path handle mapping
    path_handle_wt.push_back(0);
    path_rev_iv.push_back(0);
    path_next_id_iv.push_back(0);
    path_next_rank_iv.push_back(0);
    path_prev_id_iv.push_back(0);
    path_prev_rank_iv.push_back(0);
    // increment node count
    ++_node_count;
    // return handle
    return handle_helper::pack(new_id, 0);
}
    
/// Remove the node belonging to the given handle and all of its edges.
/// Does not update any stored paths.
/// Invalidates the destroyed handle.
/// May be called during serial for_each_handle iteration **ONLY** on the node being iterated.
/// May **NOT** be called during parallel for_each_handle iteration.
/// May **NOT** be called on the node from which edges are being followed during follow_edges.
void graph_t::destroy_handle(const handle_t& handle) {
    uint64_t offset = handle_helper::unpack_number(handle);
    id_t id = graph_id_pv.at(offset);
    // remove occs in edge lists
    // enumerate the edges
    std::vector<edge_t> edges_to_destroy;
    follow_edges(handle, false, [&edges_to_destroy,&handle,this](const handle_t& h) {
            edges_to_destroy.push_back(make_pair(handle, h)); });
    follow_edges(handle, true, [&edges_to_destroy,&handle,this](const handle_t& h) {
            edges_to_destroy.push_back(make_pair(h, handle)); });
    // and then remove them
    for (auto& edge : edges_to_destroy) {
        destroy_edge(edge);
    }
    // save the node sequence for stashing in the paths
    std::string seq = get_sequence(handle);
    // remove the sequence from seq_pv
    uint64_t seq_pv_offset = seq_bv.select1(offset);
    uint64_t length = get_length(handle);
    for (uint64_t i = seq_pv_offset; i < seq_pv_offset+length; ++i) {
        seq_pv.remove(i);
        seq_bv.remove(i);
    }
    // move the sequence of the node into each path that traverses it
    std::vector<occurrence_handle_t> occs;
    // remove reference to the node from the paths, pointing them to a new hidden node
    for_each_occurrence_on_handle(handle, [this,&occs](const occurrence_handle_t& occ) {
            occs.push_back(occ);
        });
    if (occs.size()) {
        handle_t hidden = create_hidden_handle(seq);
        for (auto& occ : occs) {
            handle_t h = get_occurrence(occ);
            if (handle_helper::unpack_bit(h)) {
                set_occurrence(occ, handle_helper::toggle_bit(hidden));
            } else {
                set_occurrence(occ, hidden);
            }
        }
    }
    // remove from path handle mapping
    do {
        path_handle_wt.remove(offset);
        path_rev_iv.remove(offset);
        path_next_id_iv.remove(offset);
        path_next_rank_iv.remove(offset);
        path_prev_id_iv.remove(offset);
        path_prev_rank_iv.remove(offset);
    } while (path_handle_wt.at(offset) != 0);
    // remove from graph_id_pv
    graph_id_pv.remove(offset);
    // from the id to handle map
    graph_id_map.erase(id);
    // and from the set of hidden nodes, if it's a member
    if (graph_id_hidden_set.count(id)) {
        graph_id_hidden_set.erase(id);
        --_hidden_count;
    }
    --_node_count;
}
    
/// Create an edge connecting the given handles in the given order and orientations.
/// Ignores existing edges.
void graph_t::create_edge(const handle_t& left, const handle_t& right) {
    //std::cerr << "create_edge from " << get_id(left) << " to " << get_id(right) << std::endl;
    //if (has_edge(left, right)) return; // do nothing if edge exists
    //std::cerr << "proceeding" << std::endl;
    handle_t left_h = left;
    handle_t right_h = right;
    if (handle_helper::unpack_bit(left_h)
        && handle_helper::unpack_bit(right_h)) {
        std::swap(left_h, right_h);
        left_h = handle_helper::toggle_bit(left_h);
        right_h = handle_helper::toggle_bit(right_h);
    }
    uint64_t left_rank = handle_helper::unpack_number(left_h);
    bool left_rev = handle_helper::unpack_bit(left_h);
    uint64_t right_rank = handle_helper::unpack_number(right_h);
    bool right_rev = handle_helper::unpack_bit(right_h);
    bool inv = (left_rev != right_rev);
    // establish the insertion value
    uint64_t right_relative = edge_to_delta(right_h, left_h);
    uint64_t left_relative = edge_to_delta(left_h, right_h);
    if (!left_rev) {
        //std::cerr << "not left rev" << std::endl;
        uint64_t edge_fwd_left_offset = edge_fwd_bv.select1(left_rank+1);
        //std::cerr << "edge fwd " << edge_fwd_left_offset << std::endl;
        edge_fwd_iv.insert(edge_fwd_left_offset, left_relative);
        edge_fwd_bv.insert(edge_fwd_left_offset, 0);
        edge_fwd_inv_bv.insert(edge_fwd_left_offset, inv);
    } else {
        //std::cerr << "left rev" << std::endl;
        uint64_t edge_rev_left_offset = edge_rev_bv.select1(left_rank+1);
        edge_rev_iv.insert(edge_rev_left_offset, left_relative);
        edge_rev_bv.insert(edge_rev_left_offset, 0);
        edge_rev_inv_bv.insert(edge_rev_left_offset, inv);
    }
    if (!right_rev) {
        //std::cerr << "not right rev" << std::endl;
        uint64_t edge_rev_right_offset = edge_rev_bv.select1(right_rank+1);
        edge_rev_iv.insert(edge_rev_right_offset, right_relative);
        edge_rev_bv.insert(edge_rev_right_offset, 0);
        edge_rev_inv_bv.insert(edge_rev_right_offset, inv);
    } else {
        //std::cerr << "right rev" << std::endl;
        uint64_t edge_fwd_right_offset = edge_fwd_bv.select1(right_rank+1);
        edge_fwd_iv.insert(edge_fwd_right_offset, right_relative);
        edge_fwd_bv.insert(edge_fwd_right_offset, 0);
        edge_fwd_inv_bv.insert(edge_fwd_right_offset, inv);
    }
    ++_edge_count;
}

uint64_t graph_t::edge_delta_to_id(uint64_t base, uint64_t delta) const {
    assert(delta != 0);
    if (delta == 1) {
        return base;
    } else if (delta % 2 == 0) {
        return base + delta/2;
    } else { //if (delta-1 % 2 == 0) {
        return base - (delta-1)/2;
    }
}

uint64_t graph_t::edge_to_delta(const handle_t& left, const handle_t& right) const {
    int64_t delta = get_id(right) - get_id(left);
    return (delta == 0 ? 1 : (delta > 0 ? 2*abs(delta) : 2*abs(delta)+1));
}

bool graph_t::has_edge(const handle_t& left, const handle_t& right) const {
    bool exists = false;
    follow_edges(left, false, [&right, &exists](const handle_t& next) {
            if (next == right) exists = true;
        });
    return exists;
}

/// Remove the edge connecting the given handles in the given order and orientations.
/// Ignores nonexistent edges.
/// Does not update any stored paths.
void graph_t::destroy_edge(const handle_t& left, const handle_t& right) {
    handle_t left_h = left;
    handle_t right_h = right;
    if (handle_helper::unpack_bit(left_h)
        && handle_helper::unpack_bit(right_h)) {
        std::swap(left_h, right_h);
        left_h = handle_helper::toggle_bit(left_h);
        right_h = handle_helper::toggle_bit(right_h);
    }
    uint64_t left_rank = handle_helper::unpack_number(left_h);
    bool left_rev = handle_helper::unpack_bit(left_h);
    uint64_t right_rank = handle_helper::unpack_number(right_h);
    bool right_rev = handle_helper::unpack_bit(right_h);
    bool inv = (left_rev != right_rev);
    // establish the insertion value
    // XXX looks strange, isn't this supposed to be ID relative?
    uint64_t right_relative = edge_to_delta(right_h, left_h);
    uint64_t left_relative = edge_to_delta(left_h, right_h);
    if (!left_rev) {
        uint64_t edge_fwd_left_offset = edge_fwd_bv.select1(left_rank);
        uint64_t edge_fwd_left_offset_erase = 0;
        for (uint64_t i = edge_fwd_left_offset+1; ; ++i) {
            uint64_t c = edge_fwd_iv.at(i);
            if (c != 0) break;
            if (c == left_relative && inv == edge_fwd_inv_bv.at(i)) {
                edge_fwd_left_offset_erase = i;
                break;
            }
        }
        if (edge_fwd_left_offset_erase) {
            edge_fwd_iv.remove(edge_fwd_left_offset_erase);
            edge_fwd_bv.remove(edge_fwd_left_offset_erase);
            edge_fwd_inv_bv.remove(edge_fwd_left_offset_erase);
        }
    } else {
        uint64_t edge_rev_left_offset = edge_rev_bv.select1(left_rank);
        uint64_t edge_rev_left_offset_erase = 0;
        for (uint64_t i = edge_rev_left_offset+1; ; ++i) {
            uint64_t c = edge_rev_iv.at(i);
            if (c != 0) break;
            if (c == left_relative && inv == edge_rev_inv_bv.at(i)) {
                edge_rev_left_offset_erase = i;
                break;
            }
        }
        if (edge_rev_left_offset_erase) {
            edge_rev_iv.remove(edge_rev_left_offset_erase);
            edge_rev_bv.remove(edge_rev_left_offset_erase);
            edge_rev_inv_bv.remove(edge_rev_left_offset_erase);
        }
    }
    if (!right_rev) {
        uint64_t edge_rev_right_offset = edge_rev_bv.select1(right_rank);
        uint64_t edge_rev_right_offset_erase = 0;
        for (uint64_t i = edge_rev_right_offset+1; ; ++i) {
            uint64_t c = edge_rev_iv.at(i);
            if (c != 0) break;
            if (c == right_relative && inv == edge_rev_inv_bv.at(i)) {
                edge_rev_right_offset_erase = i;
                break;
            }
        }
        if (edge_rev_right_offset_erase) {
            edge_rev_iv.remove(edge_rev_right_offset_erase);
            edge_rev_bv.remove(edge_rev_right_offset_erase);
            edge_rev_inv_bv.remove(edge_rev_right_offset_erase);
        }
    } else {
        uint64_t edge_fwd_right_offset = edge_fwd_bv.select1(right_rank);
        uint64_t edge_fwd_right_offset_erase = 0;
        for (uint64_t i = edge_fwd_right_offset+1; ; ++i) {
            uint64_t c = edge_fwd_iv.at(i);
            if (c != 0) break;
            if (c == right_relative && inv == edge_fwd_inv_bv.at(i)) {
                edge_fwd_right_offset_erase = i;
                break;
            }
        }
        if (edge_fwd_right_offset_erase) {
            edge_fwd_iv.remove(edge_fwd_right_offset_erase);
            edge_fwd_bv.remove(edge_fwd_right_offset_erase);
            edge_fwd_inv_bv.remove(edge_fwd_right_offset_erase);
        }
    }
    --_edge_count;
}
        
/// Remove all nodes and edges. Does not update any stored paths.
void graph_t::clear(void) {
    wt_str null_wt;
    suc_bv null_bv;
    lciv_iv null_iv;
    dyn::packed_vector null_pv;
    _max_node_id = 0;
    _min_node_id = 0;
    _node_count = 0;
    _edge_count = 0;
    _path_count = 0;
    _path_handle_next = 0;
    graph_id_pv = null_pv;
    graph_id_map.clear();
    edge_fwd_iv = null_iv;
    edge_fwd_bv = null_bv;
    edge_fwd_inv_bv = null_bv;
    edge_rev_iv = null_iv;
    edge_rev_bv = null_bv;
    edge_rev_inv_bv = null_bv;
    seq_pv = null_pv;
    seq_bv = null_bv;
    path_handle_wt = null_wt;
    path_rev_iv = null_iv;
    path_next_id_iv = null_iv;
    path_next_rank_iv = null_iv;
    path_prev_id_iv = null_iv;
    path_prev_rank_iv = null_iv;
    path_metadata_map.clear();
    path_name_map.clear();
}
    
/// Swap the nodes corresponding to the given handles, in the ordering used
/// by for_each_handle when looping over the graph. Other handles to the
/// nodes being swapped must not be invalidated. If a swap is made while
/// for_each_handle is running, it affects the order of the handles
/// traversed during the current traversal (so swapping an already seen
/// handle to a later handle's position will make the seen handle be visited
/// again and the later handle not be visited at all).
void graph_t::swap_handles(const handle_t& a, const handle_t& b) {
    assert(false);
}
    
/// Alter the node that the given handle corresponds to so the orientation
/// indicated by the handle becomes the node's local forward orientation.
/// Rewrites all edges pointing to the node and the node's sequence to
/// reflect this. Invalidates all handles to the node (including the one
/// passed). Returns a new, valid handle to the node in its new forward
/// orientation. Note that it is possible for the node's ID to change.
/// Updates all stored paths. May change the ordering of the underlying
/// graph.
handle_t graph_t::apply_orientation(const handle_t& handle) {
    // do nothing if we're already in the right orientation
    if (!handle_helper::unpack_bit(handle)) return handle;
    // store edges
    vector<handle_t> edges_fwd;
    vector<handle_t> edges_rev;
    follow_edges(handle, false, [&](const handle_t& h) {
            edges_fwd.push_back(h);
        });
    follow_edges(handle, true, [&](const handle_t& h) {
            edges_rev.push_back(h);
        });
    // save the sequence's reverse complement, which we will use to add the new handle
    const std::string seq = (handle_helper::unpack_bit(handle)
                             ? reverse_complement(get_sequence(handle))
                             : get_sequence(handle));
    // we have the technology. we can rebuild it.
    handle_t new_handle = create_handle(seq);
    // get the path set
    vector<occurrence_handle_t> occurrences;
    for_each_occurrence_on_handle(handle, [&](const occurrence_handle_t& occ) {
            occurrences.push_back(occ);
        });
    for (auto& occ : occurrences) set_occurrence(occ, new_handle);
    // destroy the old handle! (done here so as to not invalidate path occurrence iteration)
    destroy_handle(handle);
    // reconnect it to the graph
    for (auto& h : edges_fwd) {
        create_edge(handle_helper::toggle_bit(new_handle), h);
    }
    for (auto& h : edges_rev) {
        create_edge(h, handle_helper::toggle_bit(new_handle));
    }
    return new_handle;
}
    
/// Split a handle's underlying node at the given offsets in the handle's
/// orientation. Returns all of the handles to the parts. Other handles to
/// the node being split may be invalidated. The split pieces stay in the
/// same local forward orientation as the original node, but the returned
/// handles come in the order and orientation appropriate for the handle
/// passed in.
/// Updates stored paths.
std::vector<handle_t> graph_t::divide_handle(const handle_t& handle, const std::vector<size_t>& offsets) {
    // convert the offsets to the forward strand, if needed
    std::vector<size_t> fwd_offsets;
    size_t length = get_length(handle);
    if (handle_helper::unpack_bit(handle)) {
        for (auto& o : offsets) fwd_offsets.push_back(length-o);
    } else {
        for (auto& o : offsets) fwd_offsets.push_back(o);
    }
    std::sort(fwd_offsets.begin(), fwd_offsets.end());
    handle_t fwd_handle = handle_helper::unpack_bit(handle) ? handle_helper::toggle_bit(handle) : handle;
    // break it into the given pieces by building up the new node sequences
    std::string seq = get_sequence(fwd_handle);
    std::vector<std::string> seqs;
    for (uint64_t i = 0; i < offsets.size()-1; ++i) {
        seqs.push_back(seq.substr(offsets[i], offsets[i+1]-offsets[i]));
    }
    // make the handles
    std::vector<handle_t> handles;
    for (auto& s : seqs) {
        handles.push_back(create_handle(s));
    }
    // and record their reverse, for use in path fixup
    std::vector<handle_t> rev_handles;
    for (auto& h : handles) {
        rev_handles.push_back(handle_helper::toggle_bit(h));
    }
    std::reverse(rev_handles.begin(), rev_handles.end());
    // connect the pieces head to tail
    for (uint64_t i = 0; i < handles.size()-1; ++i) {
        create_edge(handles[i], handles[i+1]);
    }
    // collect the handle's path context
    vector<occurrence_handle_t> occurrences;
    for_each_occurrence_on_handle(handle, [&](const occurrence_handle_t& occ) {
            occurrences.push_back(occ);
        });
    // replace path occurrences with the new handles
    for (auto& occ : occurrences) {
        handle_t h = get_occurrence(occ);
        if (handle_helper::unpack_bit(h)) {
            replace_occurrence(occ, rev_handles);
        } else {
            replace_occurrence(occ, handles);
        }
    }
    // collect the context of the forward handle
    vector<handle_t> edges_fwd;
    vector<handle_t> edges_rev;
    follow_edges(fwd_handle, false, [&](const handle_t& h) {
            edges_fwd.push_back(h);
        });
    follow_edges(fwd_handle, true, [&](const handle_t& h) {
            edges_rev.push_back(h);
        });
    // destroy the handle
    destroy_handle(fwd_handle);
    // connect the ends to the previous context
    for (auto& h : edges_rev) create_edge(h, handles.front());
    for (auto& h : edges_rev) create_edge(handles.back(), h);
    return handles;
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
    if (get_occurrence_count(path) == 0) return; // nothing to do
    // select everything with that handle in the path_handle_wt
    dyn::packed_vector path_pv;
    for_each_occurrence_in_path(path, [this,&path_pv](const occurrence_handle_t& occ) {
            path_pv.append(as_integers(occ)[1]);
        });
    // now destroy the path
    for (uint64_t i = 0; i < path_pv.size(); ++i) {
        destroy_path_handle_records(path_pv.at(i));
    }
    path_name_map.erase(get_path_name(path));
    path_metadata_map.erase(as_integer(path));
    --_path_count;
}

void graph_t::destroy_path_handle_records(uint64_t i) {
    path_handle_wt.remove(i);
    path_rev_iv.remove(i);
    path_next_id_iv.remove(i);
    path_next_rank_iv.remove(i);
    path_prev_id_iv.remove(i);
    path_prev_rank_iv.remove(i);
}

/**
 * Create a path with the given name. The caller must ensure that no path
 * with the given name exists already, or the behavior is undefined.
 * Returns a handle to the created empty path. Handles to other paths must
 * remain valid.
 */
path_handle_t graph_t::create_path_handle(const std::string& name) {
    path_handle_t path = as_path_handle(_path_handle_next++);
    path_name_map[name] = as_integer(path);
    auto& p = path_metadata_map[as_integer(path)]; // set empty record
    occurrence_handle_t occ;
    as_integers(occ)[0] = 0;
    as_integers(occ)[1] = 0;
    p.first = occ;
    p.last = occ;
    p.length = 0;
    p.name = name;
    ++_path_count;
    return path;
}

occurrence_handle_t graph_t::create_occurrence(const path_handle_t& path, const handle_t& handle) {
    // where are we going to insert?
    uint64_t rank_on_handle = get_occurrence_count(handle);
    // build our occurrence
    occurrence_handle_t occ;
    as_integers(occ)[0] = handle_helper::unpack_number(handle);
    as_integers(occ)[1] = rank_on_handle;
    // find our insertion point
    uint64_t i = occurrence_rank(occ);
    // add reference to the path handle mapping
    path_handle_wt.insert(i, as_integer(path)+1);
    // record our handle orientation
    path_rev_iv.insert(i, handle_helper::unpack_bit(handle));
    // pad the next step
    path_next_id_iv.insert(i, path_end_marker);
    path_next_rank_iv.insert(i, 0);
    // pad the previous step
    path_prev_id_iv.insert(i, path_begin_marker);
    path_prev_rank_iv.insert(i, 0);
    return occ;
}

void graph_t::link_occurrences(const occurrence_handle_t& from, const occurrence_handle_t& to) {
    //std::cerr << "linking " << as_integers(from)[0] << "/" << as_integers(from)[1] << " and " << as_integers(to)[0] << "/" << as_integers(to)[1] << std::endl;
    path_handle_t path = get_path(from);
    assert(path == get_path(to));
    uint64_t i = occurrence_rank(from);
    // WHY CAN'T I EDIT A PACKED_VECTOR??????
    // and this wt is dog slow...
    path_next_id_iv[i] = edge_to_delta(get_occurrence(from), get_occurrence(to))+2;
    path_next_rank_iv[i] = as_integers(to)[1];
    uint64_t j = occurrence_rank(to);
    path_prev_id_iv[j] = edge_to_delta(get_occurrence(to), get_occurrence(from))+2;
    path_prev_rank_iv[j] = as_integers(from)[1];
}

void graph_t::destroy_occurrence(const occurrence_handle_t& occurrence_handle) {
    // erase reference to this occurrence
    if (has_previous_occurrence(occurrence_handle)) {
        auto occ = get_previous_occurrence(occurrence_handle);
        uint64_t i = occurrence_rank(occ);
        path_next_id_iv[i] = path_end_marker;
        path_next_rank_iv[i] = 0;
    }
    if (has_next_occurrence(occurrence_handle)) {
        auto occ = get_next_occurrence(occurrence_handle);
        uint64_t i = occurrence_rank(occ);
        path_prev_id_iv[i] = path_begin_marker;
        path_prev_rank_iv[i] = 0;
    }
    // update other records on this path on this node
    handle_t handle = get_occurrence(occurrence_handle);
    bool seen_curr = false;
    for_each_occurrence_on_handle(handle, [&](const occurrence_handle_t& occ) {
            if (occ == occurrence_handle) {
                if (seen_curr) {
                    decrement_rank(occ);
                } else {
                    seen_curr = true;
                }
            }
        });
    destroy_path_handle_records(occurrence_rank(occurrence_handle));
}

/**
 * Append a visit to a node to the given path. Returns a handle to the new
 * final occurrence on the path which is appended. Handles to prior
 * occurrences on the path, and to other paths, must remain valid.
 */
occurrence_handle_t graph_t::append_occurrence(const path_handle_t& path, const handle_t& to_append) {
    // get the last occurrence
    auto& p = path_metadata_map[as_integer(path)];
    // create the new occurrence
    occurrence_handle_t new_occ = create_occurrence(path, to_append);
    if (!p.length) {
        p.first = new_occ;
    } else {
        occurrence_handle_t last_occ = get_last_occurrence(path);
        // link it to the last step
        link_occurrences(last_occ, new_occ);
    }
    // point to the new last occ
    p.last = new_occ;
    // update our occurrence count
    ++p.length;
    return new_occ;
}

/// helper to handle the case where we remove an occurrence from a given path
/// on a node that has other occurrences from the same path, thus invalidating the
/// ranks used to refer to it
void graph_t::decrement_rank(const occurrence_handle_t& occurrence_handle) {
    // what is the actual rank of this occurrence?
    if (has_previous_occurrence(occurrence_handle)) {
        auto occ = get_previous_occurrence(occurrence_handle);
        // decrement the rank information
        uint64_t i = occurrence_rank(occ);
        uint64_t p = path_next_rank_iv.at(i);
        assert(p > 0);
        path_next_rank_iv[i] = p-1;
    }
    if (has_next_occurrence(occurrence_handle)) {
        auto occ = get_next_occurrence(occurrence_handle);
        // decrement the rank information
        uint64_t i = occurrence_rank(occ);
        uint64_t p = path_prev_rank_iv.at(i);
        assert(p > 0);
        path_prev_rank_iv[i] = p-1;
    }
}

/// reassign the given occurrence to the new handle
occurrence_handle_t graph_t::set_occurrence(const occurrence_handle_t& occurrence_handle, const handle_t& assign_to) {
    return replace_occurrence(occurrence_handle, { assign_to }).front();
}

std::vector<occurrence_handle_t>
graph_t::replace_occurrence(const occurrence_handle_t& occurrence_handle,
                            const std::vector<handle_t>& handles) {
    // verify path integrity
    const std::string prev_seq = get_sequence(get_occurrence(occurrence_handle));
    std::string new_seq;
    for (auto& handle : handles) new_seq.append(get_sequence(handle));
    assert(prev_seq == new_seq);
    // find the current occurrence
    handle_t curr_handle = get_occurrence(occurrence_handle);
    // assert that we are not invalidating the path
    assert(get_sequence(curr_handle) == get_sequence(assign_to));
    // we should not try to use this to reassign things to the same node
    for (auto& handle : handles) assert(curr_handle != handle);
    // get the context
    occurrence_handle_t prev_occ, next_occ;
    if (has_previous_occurrence(occurrence_handle)) {
        prev_occ = get_previous_occurrence(occurrence_handle);
    } else {
        as_integers(prev_occ)[0] = path_begin_marker;
        as_integers(prev_occ)[1] = 0;
    }
    if (has_next_occurrence(occurrence_handle)) {
        next_occ = get_next_occurrence(occurrence_handle);
    } else {
        as_integers(next_occ)[0] = path_end_marker;
        as_integers(next_occ)[1] = 0;
    }
    // get the path
    path_handle_t path = get_path(occurrence_handle);
    // destroy the current occurrence
    destroy_occurrence(occurrence_handle);
    // determine the new occurrences
    std::vector<occurrence_handle_t> new_occs;
    for (auto& handle : handles) {
        new_occs.push_back(create_occurrence(path, handle));
    }
    // link new occurrences
    for (uint64_t i = 0; i < new_occs.size()-1; ++i) {
        link_occurrences(new_occs[i], new_occs[i+1]);
    }
    // link to context
    if (as_integers(prev_occ)[0] != path_begin_marker) {
        link_occurrences(prev_occ, new_occs.front());
    }
    if (as_integers(next_occ)[0] != path_end_marker) {
        link_occurrences(new_occs.back(), next_occ);
    }
    return new_occs;
}

void graph_t::display(void) const {
    std::cerr << "------ graph state ------" << std::endl;

    std::cerr << "_max_node_id = " << _max_node_id << std::endl;
    std::cerr << "_min_node_id = " << _min_node_id << std::endl;

    std::cerr << "graph_id_pv" << "\t";
    for (uint64_t i = 0; i < graph_id_pv.size(); ++i) std::cerr << graph_id_pv.at(i) << " "; std::cerr << std::endl;
    /// Records edges of the 3' end on the forward strand, delimited by 0
    /// ordered by rank in graph_id_pv, defined by opposite rank+1 (handle)
    std::cerr << "edge_fwd_iv" << "\t";
    for (uint64_t i = 0; i < edge_fwd_iv.size(); ++i) std::cerr << edge_fwd_iv.at(i) << " "; std::cerr << std::endl;
    std::cerr << "edge_fwd_bv" << "\t";
    for (uint64_t i = 0; i < edge_fwd_bv.size(); ++i) std::cerr << edge_fwd_bv.at(i) << " "; std::cerr << std::endl;
    /// Marks inverting edges in edge_fwd_wt
    std::cerr << "edge_fwd_inv_bv" << "\t";
    for (uint64_t i = 0; i < edge_fwd_inv_bv.size(); ++i) std::cerr << edge_fwd_inv_bv.at(i) << " "; std::cerr << std::endl;
    /// Records edges of the 3' end on the reverse strand, delimited by 0,
    /// ordered by rank in graph_id_pv, defined by opposite rank+1 (handle)
    std::cerr << "edge_rev_iv" << "\t";
    for (uint64_t i = 0; i < edge_rev_iv.size(); ++i) std::cerr << edge_rev_iv.at(i) << " "; std::cerr << std::endl;
    std::cerr << "edge_rev_bv" << "\t";
    for (uint64_t i = 0; i < edge_rev_bv.size(); ++i) std::cerr << edge_rev_bv.at(i) << " "; std::cerr << std::endl;
    /// Marks inverting edges in edge_rev_wt
    std::cerr << "edge_rev_inv_bv" << "\t";
    for (uint64_t i = 0; i < edge_rev_inv_bv.size(); ++i) std::cerr << edge_rev_inv_bv.at(i) << " "; std::cerr << std::endl;
    /// Encodes all of the sequences of all nodes and all paths in the graph.
    /// The node sequences occur in the same order as in graph_iv;
    /// Node boundaries are given by 0s
    std::cerr << "seq_pv" << "\t\t";
    for (uint64_t i = 0; i < seq_pv.size(); ++i) std::cerr << seq_pv.at(i) << " "; std::cerr << std::endl;
    std::cerr << "seq_bv" << "\t\t";
    for (uint64_t i = 0; i < seq_bv.size(); ++i) std::cerr << seq_bv.at(i) << " "; std::cerr << std::endl;
    /// Ordered across the nodes in graph_id_pv, stores the path ids (1-based) at each
    /// segment in seq_wt, delimited by 0, one for each path occurrrence (node traversal).
    std::cerr << "path_handle_wt" << "\t";
    for (uint64_t i = 0; i < path_handle_wt.size(); ++i) std::cerr << path_handle_wt.at(i) << " "; std::cerr << std::endl;
    std::cerr << "path_rev_iv" << "\t";
    for (uint64_t i = 0; i < path_rev_iv.size(); ++i) std::cerr << path_rev_iv.at(i) << " "; std::cerr << std::endl;
    std::cerr << "path_next_id_iv" << "\t";
    for (uint64_t i = 0; i < path_next_id_iv.size(); ++i) {
        //std::cerr << "i is " << i << std::endl;
        uint64_t j = path_next_id_iv.at(i);
        if (j == path_begin_marker) std::cerr << "^";
        else if (j == path_end_marker) std::cerr << "$";
        else std::cerr << j;
        std::cerr << " ";
    } std::cerr << std::endl;
    std::cerr << "path_next_rn_wt" << "\t";
    for (uint64_t i = 0; i < path_next_rank_iv.size(); ++i) std::cerr << path_next_rank_iv.at(i) << " "; std::cerr << std::endl;
    std::cerr << "path_prev_id_iv" << "\t";
    for (uint64_t i = 0; i < path_prev_id_iv.size(); ++i) {
        uint64_t j = path_prev_id_iv.at(i);
        if (j == path_begin_marker) std::cerr << "^";
        else if (j == path_end_marker) std::cerr << "$";
        else std::cerr << j;
        std::cerr << " ";
    } std::cerr << std::endl;
    std::cerr << "path_prev_rn_wt" << "\t";
    for (uint64_t i = 0; i < path_prev_rank_iv.size(); ++i) std::cerr << path_prev_rank_iv.at(i) << " "; std::cerr << std::endl;
    std::cerr << "path_metadata" << "\t";
    for (auto& p : path_metadata_map) {
        std::cerr << p.first << ":" << p.second.name << ":"
                  << as_integers(p.second.first)[0] << "/" << as_integers(p.second.first)[1] << "->"
                  << as_integers(p.second.last)[0] << "/" << as_integers(p.second.last)[1] << " ";
    } std::cerr << std::endl;

    // not dumped...
    /// Stores path names in their internal order, delimited by '$'
    //dyn::wt_fmi path_name_fmi;
    /// Marks the beginning of each path name
    //dyn::suc_bv path_name_bv;

}

void graph_t::to_gfa(std::ostream& out) const {
    out << "H\tVN:Z:1.0" << std::endl;
    // for each node
    for_each_handle([&out,this](const handle_t& h) {
            out << "S\t" << get_id(h) << "\t" << get_sequence(h) << std::endl;
            // get the forward edges from this handle
            follow_edges(h, false, [&out, &h, this](const handle_t& a){
                    out << "L\t" << get_id(h) << "\t"
                        << (handle_helper::unpack_bit(h)?"-":"+")
                        << "\t" << get_id(a) << "\t"
                        << (handle_helper::unpack_bit(a)?"-":"+")
                        << "\t0M" << std::endl;
                });
        });
    for_each_path_handle([&out,this](const path_handle_t& p) {
            //occurrence_handle_t occ = get_first_occurrence(p);
            out << "P\t" << get_path_name(p) << "\t";
            for_each_occurrence_in_path(p, [this,&out](const occurrence_handle_t& occ) {
                    handle_t h = get_occurrence(occ);
                    out << get_id(h) << (handle_helper::unpack_bit(h)?"-":"+");
                    if (has_next_occurrence(occ)) out << ",";
                });
            out << "\t";
            for_each_occurrence_in_path(p, [this,&out](const occurrence_handle_t& occ) {
                    out << get_length(get_occurrence(occ)) << "M";
                    if (has_next_occurrence(occ)) out << ",";
                });
            out << std::endl;
        });
}

uint64_t graph_t::serialize(std::ostream& out) const {
    uint64_t written = 0;
    out.write((char*)&_max_node_id,sizeof(_max_node_id));
    written += sizeof(_max_node_id);
    out.write((char*)&_min_node_id,sizeof(_min_node_id));
    written += sizeof(_min_node_id);
    out.write((char*)&_node_count,sizeof(_node_count));
    written += sizeof(_node_count);
    out.write((char*)&_edge_count,sizeof(_edge_count));
    written += sizeof(_edge_count);
    out.write((char*)&_path_count,sizeof(_path_count));
    written += sizeof(_path_count);
    out.write((char*)&_path_handle_next,sizeof(_path_handle_next));
    written += sizeof(_path_handle_next);
    written += graph_id_pv.serialize(out);
    size_t i = graph_id_map.size();
    out.write((char*)&i,sizeof(size_t));
    written += sizeof(size_t);
    for (auto& p : graph_id_map) {
        out.write((char*)&p.first,sizeof(p.first));
        written += sizeof(p.first);
        out.write((char*)&p.second,sizeof(p.second));
        written += sizeof(p.second);
    }
    written += edge_fwd_iv.serialize(out);
    written += edge_fwd_bv.serialize(out);
    written += edge_fwd_inv_bv.serialize(out);
    written += edge_rev_iv.serialize(out);
    written += edge_rev_bv.serialize(out);
    written += edge_rev_inv_bv.serialize(out);
    written += seq_pv.serialize(out);
    written += seq_bv.serialize(out);
    written += path_handle_wt.serialize(out);
    written += path_rev_iv.serialize(out);
    written += path_next_id_iv.serialize(out);
    written += path_next_rank_iv.serialize(out);
    written += path_prev_id_iv.serialize(out);
    written += path_prev_rank_iv.serialize(out);
    i = path_metadata_map.size();
    out.write((char*)&i,sizeof(size_t));
    written += sizeof(size_t);
    for (auto& p : path_metadata_map) {
        out.write((char*)&p.first,sizeof(p.first));
        written += sizeof(p.first);
        auto& m = p.second;
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

void graph_t::load(std::istream& in) {
    //uint64_t written = 0;
    in.read((char*)&_max_node_id,sizeof(_max_node_id));
    in.read((char*)&_min_node_id,sizeof(_min_node_id));
    in.read((char*)&_node_count,sizeof(_node_count));
    in.read((char*)&_edge_count,sizeof(_edge_count));
    in.read((char*)&_path_count,sizeof(_path_count));
    in.read((char*)&_path_handle_next,sizeof(_path_handle_next));
    graph_id_pv.load(in);
    size_t i;
    in.read((char*)&i,sizeof(size_t));
    for (size_t j = 0; j < i; ++j) {
        uint64_t k, v;
        in.read((char*)&k,sizeof(uint64_t));
        in.read((char*)&v,sizeof(uint64_t));
        graph_id_map[k] = v;
    }
    edge_fwd_iv.load(in);
    edge_fwd_bv.load(in);
    edge_fwd_inv_bv.load(in);
    edge_rev_iv.load(in);
    edge_rev_bv.load(in);
    edge_rev_inv_bv.load(in);
    seq_pv.load(in);
    seq_bv.load(in);
    path_handle_wt.load(in);
    path_rev_iv.load(in);
    path_next_id_iv.load(in);
    path_next_rank_iv.load(in);
    path_prev_id_iv.load(in);
    path_prev_rank_iv.load(in);
    in.read((char*)&i,sizeof(size_t));
    for (size_t j = 0; j < i; ++j) {
        uint64_t k;
        in.read((char*)&k,sizeof(uint64_t));
        auto& m = path_metadata_map[k];
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
    for (size_t j = 0; j < i; ++j) {
        uint64_t s;
        in.read((char*)&s,sizeof(size_t));
        char k[s+1]; k[s] = '\0';
        in.read(k,s);
        uint64_t v;
        in.read((char*)&v,sizeof(uint64_t));
        path_name_map[string(k)] = v;
    }
}

}
