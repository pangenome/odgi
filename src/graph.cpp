//
//  graph.cpp
//  

#include "graph.hpp"

namespace dankgraph {

graph_t::graph_t(void) {
    // set up initial delimiters
    seq_bv.push_back(1);
    edge_fwd_wt.push_back(0);
    edge_fwd_inv_bv.push_back(0);
    edge_rev_wt.push_back(0);
    edge_rev_inv_bv.push_back(0);
    path_handle_wt.push_back(path_handle_wt_end_marker);
    path_seq_wt.push_back(0);
    path_name_fmi.extend('$');
    path_name_pv.push_back(0);
    path_name_bv.push_back(1);
}

graph_t::~graph_t(void) { }

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
        uint64_t edges_begin = edge_fwd_wt.select(offset, 0)+1;
        for (uint64_t i = edges_begin; ; ++i) {
            uint64_t x = edge_fwd_wt.at(i);
            if (x==0) break; // end of record
            uint64_t id = edge_delta_to_id(get_id(handle), x);
            bool inv = edge_fwd_inv_bv.at(i);
            handle_t handle = get_handle(id, (inv ? !is_rev : is_rev));
            result &= iteratee(handle);
            if (!result) break;
        }
    } else {
        assert(go_left && !is_rev || !go_left && is_rev);
        uint64_t edges_begin = edge_rev_wt.select(offset, 0)+1;
        for (uint64_t i = edges_begin; ; ++i) {
            uint64_t x = edge_rev_wt.at(i);
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
    std::vector<uint64_t> occs = path_name_fmi.locate(query_v);
    assert(occs.size() <= 1);
    if (occs.empty()) return false;
    uint64_t offset = occs.front();
    path_handle_t path = as_path_handle(path_name_bv.rank1(offset));
    return get_occurrence_count(path) > 0;
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
    std::string name;
    uint64_t name_begin = path_name_bv.select1(as_integer(path_handle))+1;
    uint64_t name_end = path_name_bv.select1(as_integer(path_handle)+1);
    for (uint64_t i = name_begin; i < name_end; ++i) {
        name += path_name_pv.at(i);
    }
    return name;
}
    
/// Returns the number of node occurrences in the path
size_t graph_t::get_occurrence_count(const path_handle_t& path_handle) const {
    uint64_t path_begin = path_handle_wt.select(as_integer(path_handle), path_handle_wt_end_marker)+1;
    uint64_t path_end = path_handle_wt.select(as_integer(path_handle)+1, path_handle_wt_end_marker);
    return path_end - path_begin;
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
    handle_t handle_fwd = handle;
    handle_t handle_rev = handle_helper::toggle_bit(handle);
    uint64_t occs_on_handle_fwd = path_handle_wt.rank(path_handle_wt.size(), as_integer(handle_fwd));
    uint64_t occs_on_handle_rev = path_handle_wt.rank(path_handle_wt.size(), as_integer(handle_rev));
    for (uint64_t i = 0; i < occs_on_handle_fwd; ++i) {
        uint64_t o = path_handle_wt.select(i, as_integer(handle_fwd));
        occurrence_handle_t occ;
        as_integers(occ)[0] = path_handle_wt.rank(o, path_handle_wt_end_marker)-1;
        uint64_t path_begin = path_handle_wt.select(as_integers(occ)[0], path_handle_wt_end_marker);
        as_integers(occ)[1] = i - path_handle_wt.rank(path_begin, as_integer(handle_fwd))-1;
        iteratee(occ);
    }
    for (uint64_t i = 0; i < occs_on_handle_rev; ++i) {
        uint64_t o = path_handle_wt.select(i, as_integer(handle_rev));
        occurrence_handle_t occ;
        as_integers(occ)[0] = path_handle_wt.rank(o, path_handle_wt_end_marker)-1;
        uint64_t path_begin = path_handle_wt.select(as_integers(occ)[0], path_handle_wt_end_marker);
        as_integers(occ)[1] = i - path_handle_wt.rank(path_begin, as_integer(handle_rev))-1;
        iteratee(occ);
    }
}

/// Get a node handle (node ID and orientation) from a handle to an occurrence on a path
handle_t graph_t::get_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    path_handle_t path = as_path_handle(occ_handle[0]);
    uint64_t path_begin = path_handle_wt.select(as_integer(path), path_handle_wt_end_marker)+1;
    return as_handle(path_handle_wt.at(path_begin + occ_handle[1]));
}

/// Get a handle to the first occurrence in a path.
/// The path MUST be nonempty.
occurrence_handle_t graph_t::get_first_occurrence(const path_handle_t& path_handle) const {
    assert(get_occurrence_count(path_handle) > 0);
    occurrence_handle_t result;
    as_integers(result)[0] = as_integer(path_handle);
    as_integers(result)[1] = 0;
    return result;
}
    
/// Get a handle to the last occurrence in a path
/// The path MUST be nonempty.
occurrence_handle_t graph_t::get_last_occurrence(const path_handle_t& path_handle) const {
    uint64_t path_begin = path_handle_wt.select(as_integer(path_handle), path_handle_wt_end_marker)+1;
    uint64_t path_end = path_handle_wt.select(as_integer(path_handle)+1, path_handle_wt_end_marker);
    assert(path_end > path_begin);
    occurrence_handle_t result;
    as_integers(result)[0] = as_integer(path_handle);
    as_integers(result)[1] = path_end - path_begin;
    return result;
}
    
/// Returns true if the occurrence is not the last occurence on the path, else false
bool graph_t::has_next_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    return get_occurrence_count(as_path_handle(occ_handle[0])) > occ_handle[1]+1;
}
    
/// Returns true if the occurrence is not the first occurence on the path, else false
bool graph_t::has_previous_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    return get_occurrence_count(as_path_handle(occ_handle[0])) && occ_handle[1] > 0;
}

/// Returns a handle to the next occurrence on the path, which must exist
occurrence_handle_t graph_t::get_next_occurrence(const occurrence_handle_t& occurrence_handle) const {
    occurrence_handle_t occ = occurrence_handle;
    as_integers(occ)[1]++;
    return occ;
}

/// Returns a handle to the previous occurrence on the path
occurrence_handle_t graph_t::get_previous_occurrence(const occurrence_handle_t& occurrence_handle) const {
    occurrence_handle_t occ = occurrence_handle;
    as_integers(occ)[1]--;
    return occ;
}
    
/// Returns a handle to the path that an occurrence is on
path_handle_t graph_t::get_path_handle_of_occurrence(const occurrence_handle_t& occurrence_handle) const {
    return as_path_handle(as_integers(occurrence_handle)[0]);
}
    
/// Returns the 0-based ordinal rank of a occurrence on a path
size_t graph_t::get_ordinal_rank_of_occurrence(const occurrence_handle_t& occurrence_handle) const {
    return as_integers(occurrence_handle)[1];
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
    edge_fwd_wt.push_back(0);
    edge_fwd_inv_bv.push_back(0);
    edge_rev_wt.push_back(0);
    edge_rev_inv_bv.push_back(0);
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
    // remove from graph_id_pv
    id_t id = graph_id_pv.at(offset);
    graph_id_pv.remove(offset);
    graph_id_map.erase(id);
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
    // remove reference to the node from the paths
    for_each_occurrence_on_handle(handle, [&](const occurrence_handle_t& occ) {
            unlink_occurrence(occ);
        });
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
        uint64_t edge_fwd_left_offset = edge_fwd_wt.select(left_rank+1, 0);
        //std::cerr << "edge fwd " << edge_fwd_left_offset << std::endl;
        edge_fwd_wt.insert(edge_fwd_left_offset, left_relative);
        edge_fwd_inv_bv.insert(edge_fwd_left_offset, inv);
    } else {
        //std::cerr << "left rev" << std::endl;
        uint64_t edge_rev_left_offset = edge_rev_wt.select(left_rank+1, 0);
        edge_rev_wt.insert(edge_rev_left_offset, left_relative);
        edge_rev_inv_bv.insert(edge_rev_left_offset, inv);
    }
    if (!right_rev) {
        //std::cerr << "not right rev" << std::endl;
        uint64_t edge_rev_right_offset = edge_rev_wt.select(right_rank+1, 0);
        edge_rev_wt.insert(edge_rev_right_offset, right_relative);
        edge_rev_inv_bv.insert(edge_rev_right_offset, inv);
    } else {
        //std::cerr << "right rev" << std::endl;
        uint64_t edge_fwd_right_offset = edge_fwd_wt.select(right_rank+1, 0);
        edge_fwd_wt.insert(edge_fwd_right_offset, right_relative);
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
    uint64_t right_relative = edge_to_delta(right_h, left_h);
    uint64_t left_relative = edge_to_delta(left_h, right_h);
    if (!left_rev) {
        uint64_t edge_fwd_left_offset = edge_fwd_wt.select(left_rank, 0);
        uint64_t edge_fwd_left_offset_erase = 0;
        for (uint64_t i = edge_fwd_left_offset+1; ; ++i) {
            uint64_t c = edge_fwd_wt.at(i);
            if (c != 0) break;
            if (c == left_relative && inv == edge_fwd_inv_bv.at(i)) {
                edge_fwd_left_offset_erase = i;
                break;
            }
        }
        if (edge_fwd_left_offset_erase) {
            edge_fwd_wt.remove(edge_fwd_left_offset_erase);
            edge_fwd_inv_bv.remove(edge_fwd_left_offset_erase);
        }
    } else {
        uint64_t edge_rev_left_offset = edge_rev_wt.select(left_rank, 0);
        uint64_t edge_rev_left_offset_erase = 0;
        for (uint64_t i = edge_rev_left_offset+1; ; ++i) {
            uint64_t c = edge_rev_wt.at(i);
            if (c != 0) break;
            if (c == left_relative && inv == edge_rev_inv_bv.at(i)) {
                edge_rev_left_offset_erase = i;
                break;
            }
        }
        if (edge_rev_left_offset_erase) {
            edge_rev_wt.remove(edge_rev_left_offset_erase);
            edge_rev_inv_bv.remove(edge_rev_left_offset_erase);
        }
    }
    if (!right_rev) {
        uint64_t edge_rev_right_offset = edge_rev_wt.select(right_rank, 0);
        uint64_t edge_rev_right_offset_erase = 0;
        for (uint64_t i = edge_rev_right_offset+1; ; ++i) {
            uint64_t c = edge_rev_wt.at(i);
            if (c != 0) break;
            if (c == right_relative && inv == edge_rev_inv_bv.at(i)) {
                edge_rev_right_offset_erase = i;
                break;
            }
        }
        if (edge_rev_right_offset_erase) {
            edge_rev_wt.remove(edge_rev_right_offset_erase);
            edge_rev_inv_bv.remove(edge_rev_right_offset_erase);
        }
    } else {
        uint64_t edge_fwd_right_offset = edge_fwd_wt.select(right_rank, 0);
        uint64_t edge_fwd_right_offset_erase = 0;
        for (uint64_t i = edge_fwd_right_offset+1; ; ++i) {
            uint64_t c = edge_fwd_wt.at(i);
            if (c != 0) break;
            if (c == right_relative && inv == edge_fwd_inv_bv.at(i)) {
                edge_fwd_right_offset_erase = i;
                break;
            }
        }
        if (edge_fwd_right_offset_erase) {
            edge_fwd_wt.remove(edge_fwd_right_offset_erase);
            edge_fwd_inv_bv.remove(edge_fwd_right_offset_erase);
        }
    }
    --_edge_count;
}
        
/// Remove all nodes and edges. Does not update any stored paths.
void graph_t::clear(void) {
    wt_str null_wt;
    suc_bv null_bv;
    dyn::packed_vector null_pv;
    dyn::wt_fmi null_fmi;
    _max_node_id = 0;
    _min_node_id = 0;
    graph_id_pv = null_pv;
    edge_fwd_wt = null_wt;
    edge_fwd_inv_bv = null_bv;
    edge_rev_wt = null_wt;
    edge_rev_inv_bv = null_bv;
    seq_pv = null_pv;
    seq_bv = null_bv;
    path_handle_wt = null_wt;
    path_name_fmi = null_fmi;
    path_name_bv = null_bv;
    _node_count = 0;
    _edge_count = 0;
    _path_count = 0;
    _path_handle_next = 0;
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
    uint64_t path_begin = path_handle_wt.select(as_integer(path), path_handle_wt_end_marker)+1;
    uint64_t path_end = path_handle_wt.select(as_integer(path)+1, path_handle_wt_end_marker);
    if (path_begin - path_end == 0) return; // nothing to do
    // destroy the unlinked sequences
    for_each_occurrence_in_path(path, [this,&path_begin](const occurrence_handle_t& occ) {
            if (is_unlinked(occ)) {
                uint64_t occ_offset = path_begin + as_integers(occ)[1];
                uint64_t unlinked_rank = path_handle_wt.rank(occ_offset, path_handle_wt_unlinked_marker);
                uint64_t i = path_seq_wt.select(unlinked_rank, 0)+1;
                while (path_seq_wt.at(i) != 0) {
                    path_seq_wt.remove(i);
                }
                path_seq_wt.remove(i); // remove trailing delimiter 0
            }
        });
    // now destroy the path
    for (uint64_t i = 0; i < path_end - path_begin; ++i) {
        path_handle_wt.remove(path_begin);
    }
    --_path_count;
}

/**
 * Create a path with the given name. The caller must ensure that no path
 * with the given name exists already, or the behavior is undefined.
 * Returns a handle to the created empty path. Handles to other paths must
 * remain valid.
 */
path_handle_t graph_t::create_path_handle(const std::string& name) {
    for (auto c : name) path_name_fmi.extend(c);
    path_name_fmi.extend('$');
    for (auto c : name) path_name_pv.push_back(c);
    path_name_pv.push_back(0);
    for (auto c : name) path_name_bv.push_back(0);
    path_name_bv.push_back(1);
    path_handle_wt.push_back(path_handle_wt_end_marker);
    ++_path_count;
    return as_path_handle(_path_handle_next++);
}
    
/**
 * Append a visit to a node to the given path. Returns a handle to the new
 * final occurrence on the path which is appended. Handles to prior
 * occurrences on the path, and to other paths, must remain valid.
 */
occurrence_handle_t graph_t::append_occurrence(const path_handle_t& path, const handle_t& to_append) {
    uint64_t path_begin = path_handle_wt.select(as_integer(path), path_handle_wt_end_marker)+1;
    uint64_t path_end = path_handle_wt.select(as_integer(path)+1, path_handle_wt_end_marker);
    path_handle_wt.insert(path_end, as_integer(to_append));
    occurrence_handle_t occ;
    as_integers(occ)[0] = as_integer(path);
    as_integers(occ)[1] = path_end - path_begin;
    return occ;
}

bool graph_t::is_unlinked(const occurrence_handle_t& occurrence_handle) {
    uint64_t path_begin = path_handle_wt.select(as_integers(occurrence_handle)[0], path_handle_wt_end_marker)+1;
    uint64_t occ_offset = path_begin + as_integers(occurrence_handle)[1];
    return path_handle_wt.at(occ_offset) == path_handle_wt_unlinked_marker;
}

void graph_t::unlink_occurrence(const occurrence_handle_t& occurrence_handle) {
    // TODO this could be more memory efficient in cases of many paths if we were
    // to share the virtual "node" we create among multiple unlinked path occurrences
    const std::string seq = get_sequence(get_occurrence(occurrence_handle));
    uint64_t path_begin = path_handle_wt.select(as_integers(occurrence_handle)[0], path_handle_wt_end_marker)+1;
    uint64_t occ_offset = path_begin + as_integers(occurrence_handle)[1];
    path_handle_wt.remove(occ_offset);
    path_handle_wt.insert(occ_offset, path_handle_wt_unlinked_marker);
    uint64_t unlinked_rank = path_handle_wt.rank(occ_offset, path_handle_wt_unlinked_marker);
    uint64_t i = path_seq_wt.select(unlinked_rank, 0)+1;
    path_seq_wt.insert(i, 0); // insert the 0
    // write in reverse
    for (uint64_t j = seq.size()-1; j >= 0; --j) {
        path_seq_wt.insert(i, dna_as_int(seq.at(j)));
    }
}

std::string graph_t::get_occurrence_sequence(const occurrence_handle_t& occurrence_handle) {
    uint64_t path_begin = path_handle_wt.select(as_integers(occurrence_handle)[0], path_handle_wt_end_marker)+1;
    uint64_t occ_offset = path_begin + as_integers(occurrence_handle)[1];
    uint64_t val = path_handle_wt.at(occ_offset);
    if (val != path_handle_wt_unlinked_marker) return get_sequence(as_handle(val));
    //if (path_handle_wt.at(occ_offset) != path_handle_wt_unlinked_marker) return seq;
    uint64_t unlinked_rank = path_handle_wt.rank(occ_offset, path_handle_wt_unlinked_marker);
    uint64_t i = path_seq_wt.select(unlinked_rank, 0)+1;
    uint64_t j = path_seq_wt.select(unlinked_rank+1, 0);
    std::string seq;
    for ( ; i < j; ++i) {
        seq += int_as_dna(path_seq_wt.at(i));
    }
    return seq;
}

void graph_t::set_occurrence(const occurrence_handle_t& occurrence_handle, const handle_t& handle) {
    assert(get_occurrence_sequence(occurrence_handle) == get_sequence(handle));
    uint64_t path_begin = path_handle_wt.select(as_integers(occurrence_handle)[0], path_handle_wt_end_marker)+1;
    uint64_t occ_offset = path_begin + as_integers(occurrence_handle)[1];
    // remove unlinked sequence
    if (is_unlinked(occurrence_handle)) {
        uint64_t unlinked_rank = path_handle_wt.rank(occ_offset, path_handle_wt_unlinked_marker);
        uint64_t i = path_seq_wt.select(unlinked_rank, 0)+1;
        while (path_seq_wt.at(i) != 0) {
            path_seq_wt.remove(i);
        }
        path_seq_wt.remove(i); // remove trailing delimiter 0
    }
    path_handle_wt.remove(occ_offset);
    path_handle_wt.insert(occ_offset, as_integer(handle));
}

void graph_t::replace_occurrence(const occurrence_handle_t& occurrence_handle, const std::vector<handle_t>& handles) {
    // verify path integrity
    const std::string prev_seq = get_sequence(get_occurrence(occurrence_handle));
    std::string new_seq;
    for (auto& handle : handles) new_seq.append(get_sequence(handle));
    assert(prev_seq == new_seq);
    uint64_t path_begin = path_handle_wt.select(as_integers(occurrence_handle)[0], path_handle_wt_end_marker)+1;
    uint64_t occ_offset = path_begin + as_integers(occurrence_handle)[1];
    // remove any unlinked sequence
    if (is_unlinked(occurrence_handle)) {
        uint64_t unlinked_rank = path_handle_wt.rank(occ_offset, path_handle_wt_unlinked_marker);
        uint64_t i = path_seq_wt.select(unlinked_rank, 0)+1;
        while (path_seq_wt.at(i) != 0) {
            path_seq_wt.remove(i);
        }
        path_seq_wt.remove(i); // remove trailing delimiter 0
    }
    // make the replacement
    path_handle_wt.remove(occ_offset);
    for (uint64_t i = handles.size()-1; i >= 0; --i) {
        path_handle_wt.insert(occ_offset, as_integer(handles.at(i)));
    }
}

void graph_t::display(void) const {
    std::cerr << "------ graph state ------" << std::endl;

    std::cerr << "_max_node_id = " << _max_node_id << std::endl;
    std::cerr << "_min_node_id = " << _min_node_id << std::endl;

    std::cerr << "graph_id_pv" << "\t";
    for (uint64_t i = 0; i < graph_id_pv.size(); ++i) std::cerr << graph_id_pv.at(i) << " "; std::cerr << std::endl;
    /// Records edges of the 3' end on the forward strand, delimited by 0
    /// ordered by rank in graph_id_pv, defined by opposite rank+1 (handle)
    std::cerr << "edge_fwd_wt" << "\t";
    for (uint64_t i = 0; i < edge_fwd_wt.size(); ++i) std::cerr << edge_fwd_wt.at(i) << " "; std::cerr << std::endl;
    /// Marks inverting edges in edge_fwd_wt
    std::cerr << "edge_fwd_inv_bv" << "\t";
    for (uint64_t i = 0; i < edge_fwd_inv_bv.size(); ++i) std::cerr << edge_fwd_inv_bv.at(i) << " "; std::cerr << std::endl;
    /// Records edges of the 3' end on the reverse strand, delimited by 0,
    /// ordered by rank in graph_id_pv, defined by opposite rank+1 (handle)
    std::cerr << "edge_rev_wt" << "\t";
    for (uint64_t i = 0; i < edge_rev_wt.size(); ++i) std::cerr << edge_rev_wt.at(i) << " "; std::cerr << std::endl;
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
    for (uint64_t i = 0; i < path_handle_wt.size(); ++i) {
        uint64_t j = path_handle_wt.at(i);
        if (j == path_handle_wt_end_marker) {
            std::cerr << "$" << " ";
        } else if (j == path_handle_wt_unlinked_marker) {
            std::cerr << "u" << " ";
        } else {
            std::cerr << j << " ";
        }
    }
    std::cerr << std::endl;
    std::cerr << "path_seq_wt" << "\t";
    for (uint64_t i = 0; i < path_seq_wt.size(); ++i) std::cerr << path_seq_wt.at(i) << " "; std::cerr << std::endl;
    std::cerr << "path_name_pv" << "\t";
    for (uint64_t i = 0; i < path_name_pv.size(); ++i) std::cerr << path_name_pv.at(i) << " "; std::cerr << std::endl;
    std::cerr << "path_name_bv" << "\t";
    for (uint64_t i = 0; i < path_name_bv.size(); ++i) std::cerr << path_name_bv.at(i) << " "; std::cerr << std::endl;

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

}
