//
//  graph.cpp
//  

#include "graph.hpp"

namespace odgi {

/// Method to check if a node exists by ID
bool graph_t::has_node(nid_t node_id) const {
    uint64_t rank = get_node_rank(node_id);
    return (rank >= node_v.size() ? false : !deleted_node_bv.at(rank));
}

/// Look up the handle for the node with the given ID in the given orientation
handle_t graph_t::get_handle(const nid_t& node_id, bool is_reverse) const {
    return number_bool_packing::pack(node_id-1, is_reverse);
}

/// Get the ID from a handle
nid_t graph_t::get_id(const handle_t& handle) const {
    return number_bool_packing::unpack_number(handle)+1;
}

/// get the backing node for a given node id
uint64_t graph_t::get_node_rank(const nid_t& node_id) const {
    return node_id - 1;
}
    
/// Get the orientation of a handle
bool graph_t::get_is_reverse(const handle_t& handle) const {
    return number_bool_packing::unpack_bit(handle);
}
    
/// Invert the orientation of a handle (potentially without getting its ID)
handle_t graph_t::flip(const handle_t& handle) const {
    return number_bool_packing::toggle_bit(handle);
}
    
/// Get the length of a node
size_t graph_t::get_length(const handle_t& handle) const {
    return node_v.at(number_bool_packing::unpack_number(handle)).sequence_size();
}

/// Get the sequence of a node, presented in the handle's local forward orientation.
std::string graph_t::get_sequence(const handle_t& handle) const {
    auto& seq = node_v.at(number_bool_packing::unpack_number(handle)).sequence();
    return (get_is_reverse(handle) ? reverse_complement(seq) : seq);
}
    
/// Loop over all the handles to next/previous (right/left) nodes. Passes
/// them to a callback which returns false to stop iterating and true to
/// continue. Returns true if we finished and false if we stopped early.
bool graph_t::follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const {
    const node_t& node = node_v.at(number_bool_packing::unpack_number(handle));
    bool is_rev = get_is_reverse(handle);
    bool result = true;
    nid_t node_id = get_id(handle);
    const std::vector<uint64_t> node_edges = node.edges();
    if (node_edges.size() == 0) return result;
    for (uint64_t i = 0; i+1 < node_edges.size(); i+=2) {
        // unpack the edge
        uint64_t other_id = edge_delta_to_id(node_id, node_edges.at(i));
        uint8_t packed_edge = node_edges.at(i+1);
        bool on_rev = edge_helper::unpack_on_rev(packed_edge);
        bool other_rev = edge_helper::unpack_other_rev(packed_edge);
        bool to_curr = edge_helper::unpack_to_curr(packed_edge);
        if (is_rev != on_rev) {
            other_rev ^= 1;
            to_curr ^= 1;
        }
        if (!go_left && !to_curr) {
            result &= iteratee(get_handle(other_id, other_rev));
        } else if (go_left && to_curr) {
            result &= iteratee(get_handle(other_id, other_rev));
        }
        if (!result) break;
    }
    return result;
}

    
/// Loop over all the nodes in the graph in their local forward
/// orientations, in their internal stored order. Stop if the iteratee
/// returns false. Can be told to run in parallel, in which case stopping
/// after a false return value is on a best-effort basis and iteration
/// order is not defined.
bool graph_t::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel) const {
    if (parallel) {
        volatile bool flag=true;
#pragma omp parallel for
        for (uint64_t i = 0; i < node_v.size(); ++i) {
            if (deleted_node_bv.at(i) == 1) continue;
            if (!flag) continue;
            bool result = iteratee(number_bool_packing::pack(i,false));
#pragma omp atomic
            flag &= result;
        }
        return flag;
    } else {
        for (uint64_t i = 0; i < node_v.size(); ++i) {
            if (deleted_node_bv.at(i) == 1) continue;
            if (!iteratee(number_bool_packing::pack(i,false))) return false;
        }
        return true;
    }
}

/// Return the number of nodes in the graph
/// TODO: can't be node_count because XG has a field named node_count.
size_t graph_t::node_size(void) const {
    return node_v.size()-_deleted_node_count;
}
    
/// Return the smallest ID in the graph, or some smaller number if the
/// smallest ID is unavailable. Return value is unspecified if the graph is empty.
nid_t graph_t::min_node_id(void) const {
    return _min_node_id;
}
    
/// Return the largest ID in the graph, or some larger number if the
/// largest ID is unavailable. Return value is unspecified if the graph is empty.
nid_t graph_t::max_node_id(void) const {
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
    size_t degree = 0;
    follow_edges(handle, go_left, [&degree](const handle_t& h) { ++degree; });
    return degree;
}
    
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
    return path_metadata_v.at(as_integer(path_handle)).name;
}
    
/// Returns the number of node occurrences in the path
size_t graph_t::get_occurrence_count(const path_handle_t& path_handle) const {
    return path_metadata_v.at(as_integer(path_handle)).length;
}

/// Returns the number of paths stored in the graph
size_t graph_t::get_path_count(void) const {
    return _path_count;
}
    
/// Execute a function on each path in the graph
bool graph_t::for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const {
    bool flag = true;
    for (uint64_t i = 0; i < _path_handle_next && flag; ++i) {
        path_handle_t path = as_path_handle(i);
        if (get_occurrence_count(path) > 0) {
            flag &= iteratee(path);
        }
    }
    return flag;
}

bool graph_t::for_each_occurrence_on_handle_impl(const handle_t& handle, const std::function<bool(const occurrence_handle_t&)>& iteratee) const {
    const node_t& node = node_v.at(number_bool_packing::unpack_number(handle));
    bool flag = true;
    nid_t handle_id = get_id(handle);
    uint64_t path_count = node.path_count();
    for (uint64_t i = 0; i < path_count; ++i) {
        occurrence_handle_t occ;
        as_integers(occ)[0] = handle_id;
        as_integers(occ)[1] = i;
        flag &= iteratee(occ);
    }
    return flag;
}

/// Returns a vector of all occurrences of a node on paths. Optionally restricts to
/// occurrences that match the handle in orientation.
std::vector<occurrence_handle_t> graph_t::occurrences_of_handle(const handle_t& handle,
                                                                bool match_orientation) const {
    std::vector<occurrence_handle_t> res;
    for_each_occurrence_on_handle(handle, [&](const occurrence_handle_t& occ) {
            handle_t h = get_occurrence(occ);
            if (!match_orientation || get_is_reverse(h) == get_is_reverse(handle)) {
                res.push_back(occ);
            }
        });
    return res;
}

size_t graph_t::get_occurrence_count(const handle_t& handle) const {
    const node_t& node = node_v.at(number_bool_packing::unpack_number(handle));
    return node.path_count();
}

/// Get a node handle (node ID and orientation) from a handle to an occurrence on a path
handle_t graph_t::get_occurrence(const occurrence_handle_t& occurrence_handle) const {
    return get_handle(as_integers(occurrence_handle)[0],
                      node_v.at(get_node_rank(as_integers(occurrence_handle)[0])).get_path_step(as_integers(occurrence_handle)[1]).is_rev());
}

/// Get a path handle (path ID) from a handle to an occurrence on a path
path_handle_t graph_t::get_path(const occurrence_handle_t& occurrence_handle) const {
    const node_t& node = node_v.at(get_node_rank(as_integers(occurrence_handle)[0]));
    return as_path_handle(node.get_path_step(as_integers(occurrence_handle)[1]).path_id());
}

/// Get a handle to the first occurrence in a path.
/// The path MUST be nonempty.
occurrence_handle_t graph_t::get_first_occurrence(const path_handle_t& path_handle) const {
    return path_metadata_v.at(as_integer(path_handle)).first;
}
    
/// Get a handle to the last occurrence in a path
/// The path MUST be nonempty.
occurrence_handle_t graph_t::get_last_occurrence(const path_handle_t& path_handle) const {
    return path_metadata_v.at(as_integer(path_handle)).last;
}
    
/// Returns true if the occurrence is not the last occurence on the path, else false
bool graph_t::has_next_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const node_t& node = node_v.at(get_node_rank(as_integers(occurrence_handle)[0]));
    return node.get_path_step(as_integers(occurrence_handle)[1]).next_id() != path_end_marker;
}
    
/// Returns true if the occurrence is not the first occurence on the path, else false
bool graph_t::has_previous_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const node_t& node = node_v.at(get_node_rank(as_integers(occurrence_handle)[0]));
    return node.get_path_step(as_integers(occurrence_handle)[1]).prev_id() != path_begin_marker;
}

/// Returns a handle to the next occurrence on the path, which must exist
occurrence_handle_t graph_t::get_next_occurrence(const occurrence_handle_t& occurrence_handle) const {
    nid_t curr_id = as_integers(occurrence_handle)[0];
    occurrence_handle_t occ;
    const node_t& node = node_v.at(get_node_rank(as_integers(occurrence_handle)[0]));
    auto& step = node.get_path_step(as_integers(occurrence_handle)[1]);
    as_integers(occ)[0] = edge_delta_to_id(curr_id, step.next_id()-2);
    as_integers(occ)[1] = step.next_rank();
    return occ;
}

/// Returns a handle to the previous occurrence on the path
occurrence_handle_t graph_t::get_previous_occurrence(const occurrence_handle_t& occurrence_handle) const {
    nid_t curr_id = as_integers(occurrence_handle)[0];
    occurrence_handle_t occ;
    const node_t& node = node_v.at(get_node_rank(as_integers(occurrence_handle)[0]));
    auto& step = node.get_path_step(as_integers(occurrence_handle)[1]);
    as_integers(occ)[0] = edge_delta_to_id(curr_id, step.prev_id()-2);
    as_integers(occ)[1] = step.prev_rank();
    return occ;
}

path_handle_t graph_t::get_path_handle_of_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const node_t& node = node_v.at(get_node_rank(as_integers(occurrence_handle)[0]));
    auto& step = node.get_path_step(as_integers(occurrence_handle)[1]);
    return as_path_handle(step.path_id());
}
    
////////////////////////////////////////////////////////////////////////////
// Additional optional interface with a default implementation
////////////////////////////////////////////////////////////////////////////

/// Returns true if the given path is empty, and false otherwise
bool graph_t::is_empty(const path_handle_t& path_handle) const {
    return get_occurrence_count(path_handle) == 0;
}

/**
 * This is the interface for a handle graph that supports modification.
 */
/*
 * Note: All operations may invalidate path handles and occurrence handles.
 */

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
    
/// Create a new node with the given sequence and return the handle.
handle_t graph_t::create_handle(const std::string& sequence) {
    // get first deleted node to recycle
    if (_deleted_node_count) {
        return create_handle(sequence, deleted_node_bv.select1(0)+1);
    } else {
        return create_handle(sequence, node_v.size()+1);
    }
}

handle_t graph_t::create_hidden_handle(const std::string& sequence) {
    // get node id as max+1
    handle_t handle = create_handle(sequence);
    nid_t id = get_id(handle);
    graph_id_hidden_set.insert(id);
    ++_hidden_count;
    return handle;
}

/// Create a new node with the given id and sequence, then return the handle.
handle_t graph_t::create_handle(const std::string& sequence, const nid_t& id) {
    assert(sequence.size());
    assert(id > 0);
    if (id > node_v.size()) {
        uint64_t to_add = id - node_v.size(); // + 1e5;
        uint64_t old_size = node_v.size();
        // realloc
        node_v.resize((uint64_t)id);
        // mark empty nodes
        for (uint64_t i = 0; i < to_add; ++i) {
            // insert before final delimiter
            deleted_node_bv.insert(old_size, 1);
            ++_deleted_node_count;
        }
    }
    // update min/max node ids
    _max_node_id = max(id, _max_node_id);
    if (_min_node_id) {
        _min_node_id = (uint64_t)min(id, _min_node_id);
    } else {
        _min_node_id = id;
    }
    // add to node vector
    uint64_t handle_rank = (uint64_t)id-1;
    // set its values
    auto& node = node_v[handle_rank];
    node.set_sequence(sequence);
    // it's not deleted
    deleted_node_bv[handle_rank] = 0;
    --_deleted_node_count;
    // increment node count
    ++_node_count;
    // return handle
    return number_bool_packing::pack(handle_rank, 0);
}
    
/// Remove the node belonging to the given handle and all of its edges.
/// Does not update any stored paths.
/// Invalidates the destroyed handle.
/// May be called during serial for_each_handle iteration **ONLY** on the node being iterated.
/// May **NOT** be called during parallel for_each_handle iteration.
/// May **NOT** be called on the node from which edges are being followed during follow_edges.
void graph_t::destroy_handle(const handle_t& handle) {
    handle_t fwd_handle = get_is_reverse(handle) ? flip(handle) : handle;
    uint64_t id = get_id(handle);
    // remove occs in edge lists
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
    // save the node sequence for stashing in the paths
    std::string seq = get_sequence(handle);
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
            if (get_is_reverse(h)) {
                set_occurrence(occ, flip(hidden));
            } else {
                set_occurrence(occ, hidden);
            }
        }
    }
    // remove from the graph by hiding it (compaction later)
    auto& node = node_v[number_bool_packing::unpack_number(handle)];
    node.clear();
    deleted_node_bv[number_bool_packing::unpack_number(handle)] = 1;
    // and from the set of hidden nodes, if it's a member
    if (graph_id_hidden_set.count(id)) {
        graph_id_hidden_set.erase(id);
        --_hidden_count;
    }
    --_node_count;
    ++_deleted_node_count;
    // check if we should compact our deleted nodes storage
}

/*
void graph_t::rebuild_id_handle_mapping(void) {
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
void graph_t::create_edge(const handle_t& left, const handle_t& right) {
    //if (has_edge(left, right)) return; // do nothing if edge exists
    handle_t left_h = left, right_h = right;
    canonicalize_edge(left_h, right_h);
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

    if (left_rank == right_rank) return;

    auto& right_node = node_v.at(right_rank);
    right_node.add_edge(right_relative,
                        edge_helper::pack(get_is_reverse(right_h),
                                          get_is_reverse(left_h),
                                          true));

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
    auto& node_l = node_v[number_bool_packing::unpack_number(left)];
    auto& node_r = node_v[number_bool_packing::unpack_number(right)];
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
void graph_t::destroy_edge(const handle_t& left, const handle_t& right) {
    //if (!has_edge(left, right)) return;
    //std::cerr << "remove edge " << get_id(left) << " " << get_id(right) << std::endl;
    handle_t left_h = left;
    handle_t right_h = right;
    canonicalize_edge(left_h, right_h);
    
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
    deleted_node_bv = null_bv;
    node_v.clear();
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
void graph_t::swap_handles(const handle_t& a, const handle_t& b) {
    //assert(false);
}

/// Reorder the graph's internal structure to match that given.
/// Optionally compact the id space of the graph to match the ordering, from 1->|ordering|.
void graph_t::apply_ordering(const std::vector<handle_t>& order, bool compact_ids) {
    graph_t ordered;
    // nodes
    hash_map<nid_t, nid_t> ids;
    ids.reserve(order.size());
    // establish id mapping
    if (compact_ids) {
        for (uint64_t i = 0; i < order.size(); ++i) {
            ids[get_id(order.at(i))] = i+1;
        }
    } else {
        for (uint64_t i = 0; i < order.size(); ++i) {
            auto& handle = order.at(i);
            ids[get_id(handle)] = get_id(handle);
        }
    }
    // nodes
    for (auto& handle : order) {
        ordered.create_handle(get_sequence(handle), ids[get_id(handle)]);
    }
    // edges
    for (auto& handle : order) {
        follow_edges(handle, false, [&](const handle_t& h) {
                ordered.create_edge(ordered.get_handle(ids[get_id(handle)], get_is_reverse(handle)),
                                    ordered.get_handle(ids[get_id(h)], get_is_reverse(h)));
            });
        follow_edges(flip(handle), false, [&](const handle_t& h) {
                ordered.create_edge(ordered.get_handle(ids[get_id(handle)], get_is_reverse(flip(handle))),
                                    ordered.get_handle(ids[get_id(h)], get_is_reverse(h)));
            });
    }
    // paths
    //std::cerr << "paths" << std::endl;
    for_each_path_handle([&](const path_handle_t& old_path) {
            //occurrence_handle_t occ = get_first_occurrence(p);
            path_handle_t new_path = ordered.create_path_handle(get_path_name(old_path));
            //std::cerr << get_path_name(old_path) << std::endl;
            for_each_occurrence_in_path(old_path, [&](const occurrence_handle_t& occ) {
                    handle_t old_handle = get_occurrence(occ);
                    handle_t new_handle = ordered.get_handle(ids[get_id(old_handle)], get_is_reverse(old_handle));
                    ordered.append_occurrence(new_path, new_handle);
                });
        });
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
handle_t graph_t::apply_orientation(const handle_t& handle) {
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
    
    node.set_sequence(get_sequence(handle));
    node.flip_paths(path_begin_marker, path_end_marker);

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

void graph_t::set_handle_sequence(const handle_t& handle, const std::string& seq) {
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
std::vector<handle_t> graph_t::divide_handle(const handle_t& handle, const std::vector<size_t>& offsets) {
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
    vector<occurrence_handle_t> occurrences;
    for_each_occurrence_on_handle(handle, [&](const occurrence_handle_t& occ) {
            occurrences.push_back(occ);
        });
    // replace path occurrences with the new handles
    for (auto& occ : occurrences) {
        handle_t h = get_occurrence(occ);
        if (get_is_reverse(h)) {
            replace_occurrence(occ, rev_handles);
        } else {
            replace_occurrence(occ, handles);
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
    std::vector<occurrence_handle_t> path_v;
    for_each_occurrence_in_path(path, [this,&path_v](const occurrence_handle_t& occ) {
            path_v.push_back(occ);
        });
    for (auto& occ : path_v) {
        destroy_occurrence(occ);
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
    path_handle_t path = as_path_handle(_path_handle_next++);
    path_name_map[name] = as_integer(path);
    path_metadata_v.emplace_back();
    auto& p = path_metadata_v.back(); //[as_integer(path)]; // set empty record
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
    as_integers(occ)[0] = get_id(handle);
    as_integers(occ)[1] = rank_on_handle;
    auto& node = node_v.at(number_bool_packing::unpack_number(handle));
    node.add_path_step(as_integer(path), get_is_reverse(handle),
                       path_begin_marker, 0, path_end_marker, 0);
    return occ;
}

void graph_t::link_occurrences(const occurrence_handle_t& from, const occurrence_handle_t& to) {
    path_handle_t path = get_path(from);
    assert(path == get_path(to));
    const handle_t& from_handle = get_occurrence(from);
    const handle_t& to_handle = get_occurrence(to);
    const uint64_t& from_rank = as_integers(from)[1];
    const uint64_t& to_rank = as_integers(to)[1];
    node_t& from_node = node_v.at(get_node_rank(as_integers(from)[0]));
    node_t::step_t from_step = from_node.get_path_step(from_rank);
    from_step.set_next_id(edge_to_delta(from_handle, to_handle)+2);
    from_step.set_next_rank(to_rank);
    from_node.set_path_step(from_rank, from_step);
    node_t& to_node = node_v.at(get_node_rank(as_integers(to)[0]));
    node_t::step_t to_step = to_node.get_path_step(to_rank);
    to_step.set_prev_id(edge_to_delta(to_handle, from_handle)+2);
    to_step.set_prev_rank(from_rank);
    to_node.set_path_step(to_rank, to_step);
}

void graph_t::destroy_occurrence(const occurrence_handle_t& occurrence_handle) {
    // erase reference to this occurrence
    bool has_prev = has_previous_occurrence(occurrence_handle);
    bool has_next = has_next_occurrence(occurrence_handle);
    if (!has_prev && !has_next) {
        // we're about to erase the path, so we need to clean up the path metadata record
        path_handle_t path = get_path_handle_of_occurrence(occurrence_handle);
        path_name_map.erase(get_path_name(path));
        path_metadata_v[as_integer(path)] = path_metadata_t();
    } else {
        if (has_prev) {
            auto occ = get_previous_occurrence(occurrence_handle);
            node_t& occ_node = node_v.at(get_node_rank(as_integers(occ)[0]));
            uint64_t occ_rank = as_integers(occ)[1];
            node_t::step_t occ_step = occ_node.get_path_step(occ_rank);
            //std::cerr << "destroy prev links " << occ_node.id() << std::endl;
            occ_step.set_next_id(path_end_marker);
            occ_step.set_next_rank(0);
            occ_node.set_path_step(occ_rank, occ_step);
        } else if (has_next) {
            auto occ = get_next_occurrence(occurrence_handle);
            auto& p = path_metadata_v[as_integer(get_path_handle_of_occurrence(occ))];
            p.first = occ;
        }
        if (has_next) {
            auto occ = get_next_occurrence(occurrence_handle);
            node_t& occ_node = node_v.at(get_node_rank(as_integers(occ)[0]));
            uint64_t occ_rank = as_integers(occ)[1];
            node_t::step_t occ_step = occ_node.get_path_step(occ_rank);
            //std::cerr << "destroy next links " << occ_node.id() << std::endl;
            occ_step.set_prev_id(path_begin_marker);
            occ_step.set_prev_rank(0);
            occ_node.set_path_step(occ_rank, occ_step);
        } else if (has_prev) {
            auto occ = get_previous_occurrence(occurrence_handle);
            auto& p = path_metadata_v[as_integer(get_path_handle_of_occurrence(occ))];
            p.last = occ;
        }
    }
    // update other records on this path on this node
    handle_t handle = get_occurrence(occurrence_handle);
    bool seen_curr = false;
    for_each_occurrence_on_handle(handle, [&](const occurrence_handle_t& occ) {
            if (seen_curr) {
                decrement_rank(occ);
            }
            if (occ == occurrence_handle) {
                seen_curr = true;
            }
        });
    node_t& curr_node = node_v.at(get_node_rank(as_integers(occurrence_handle)[0]));
    curr_node.remove_path_step(as_integers(occurrence_handle)[1]);
}

/**
 * Append a visit to a node to the given path. Returns a handle to the new
 * final occurrence on the path which is appended. Handles to prior
 * occurrences on the path, and to other paths, must remain valid.
 */
occurrence_handle_t graph_t::append_occurrence(const path_handle_t& path, const handle_t& to_append) {
    // get the last occurrence
    auto& p = path_metadata_v[as_integer(path)];
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
    //std::cerr << "in decrement rank " << as_integers(occurrence_handle)[0] << std::endl;
    if (has_previous_occurrence(occurrence_handle)) {
        auto occ = get_previous_occurrence(occurrence_handle);
        // decrement the rank information
        node_t& occ_node = node_v.at(get_node_rank(as_integers(occ)[0]));
        uint64_t occ_rank = as_integers(occ)[1];
        node_t::step_t occ_step = occ_node.get_path_step(occ_rank);
        occ_step.set_next_rank(occ_step.next_rank()-1);
        occ_node.set_path_step(occ_rank, occ_step);
    } else {
        // update path metadata
        auto& p = path_metadata_v[as_integer(get_path(occurrence_handle))];
        --as_integers(p.first)[1];
    }
    if (has_next_occurrence(occurrence_handle)) {
        auto occ = get_next_occurrence(occurrence_handle);
        node_t& occ_node = node_v.at(get_node_rank(as_integers(occ)[0]));
        uint64_t occ_rank = as_integers(occ)[1];
        node_t::step_t occ_step = occ_node.get_path_step(occ_rank);
        occ_step.set_prev_rank(occ_step.prev_rank()-1);
        occ_node.set_path_step(occ_rank, occ_step);
    } else {
        // update path metadata
        auto& p = path_metadata_v[as_integer(get_path(occurrence_handle))];
        --as_integers(p.last)[1];
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
        as_integers(prev_occ)[0] = 0;
        as_integers(prev_occ)[1] = 0;
    }
    if (has_next_occurrence(occurrence_handle)) {
        next_occ = get_next_occurrence(occurrence_handle);
    } else {
        as_integers(next_occ)[0] = 0;
        as_integers(next_occ)[1] = 0;
    }
    // get the path
    path_handle_t path = get_path(occurrence_handle);
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
    if (as_integers(prev_occ)[0]) {
        link_occurrences(prev_occ, new_occs.front());
    }
    if (as_integers(next_occ)[0]) {
        link_occurrences(new_occs.back(), next_occ);
    }
    return new_occs;
}

void graph_t::display(void) const {
    std::cerr << "------ graph state ------" << std::endl;

    std::cerr << "_max_node_id = " << _max_node_id << std::endl;
    std::cerr << "_min_node_id = " << _min_node_id << std::endl;

    //std::cerr << "graph_id_map" << "\t";
    //for (auto& k : graph_id_map) std::cerr << k.first << "->" << k.second << " "; std::cerr << std::endl;
    std::cerr << "node_v" << "\t";
    for (uint64_t i = 0; i < node_v.size(); ++i) {
        auto& node = node_v.at(i);
        nid_t node_id = i+1;
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
                      << (step.prev_id() == 0 ? "#" : std::to_string(edge_delta_to_id(node_id, step.prev_id()-2))) << ":"
                      << step.prev_rank() << ":"
                      << (step.next_id() == 1 ? "$" : std::to_string(edge_delta_to_id(node_id, step.next_id()-2))) << ":"
                      << step.next_rank() << " ";
        }
        std::cerr << " | ";
    }
    std::cerr << std::endl;
    std::cerr << "deleted_node_bv" << "\t";
    for (uint64_t i = 0; i < deleted_node_bv.size(); ++i) std::cerr << deleted_node_bv.at(i) << " "; std::cerr << std::endl;
    /// Ordered across the nodes in graph_id_iv, stores the path ids (1-based) at each
    /// segment in seq_wt, delimited by 0, one for each path occurrrence (node traversal).
    std::cerr << "path_metadata" << "\t";
    for (uint64_t q = 0; q < path_metadata_v.size(); ++q) {
        auto& p = path_metadata_v.at(q);
        std::cerr << q << ":" << p.name << ":"
                  << as_integers(p.first)[0] << "/" << as_integers(p.first)[1] << "->"
                  << as_integers(p.last)[0] << "/" << as_integers(p.last)[1] << " ";
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
                    if (as_integer(h) <= as_integer(a)) {
                        out << "L\t" << get_id(h) << "\t"
                            << (get_is_reverse(h)?"-":"+")
                            << "\t" << get_id(a) << "\t"
                            << (get_is_reverse(a)?"-":"+")
                            << "\t0M" << std::endl;
                    }
                });
            follow_edges(flip(h), false, [&out, &h, this](const handle_t& a){
                    if (as_integer(h) <= as_integer(a)) {
                        out << "L\t" << get_id(h) << "\t"
                            << (get_is_reverse(h)?"+":"-")
                            << "\t" << get_id(a) << "\t"
                            << (get_is_reverse(a)?"-":"+")
                            << "\t0M" << std::endl;
                    }
                });
        });
    for_each_path_handle([&out,this](const path_handle_t& p) {
            //occurrence_handle_t occ = get_first_occurrence(p);
            out << "P\t" << get_path_name(p) << "\t";
            for_each_occurrence_in_path(p, [this,&out](const occurrence_handle_t& occ) {
                    handle_t h = get_occurrence(occ);
                    out << get_id(h) << (get_is_reverse(h)?"-":"+");
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

uint64_t graph_t::serialize(std::ostream& out) {
    //rebuild_id_handle_mapping();
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
    out.write((char*)&_deleted_node_count,sizeof(_deleted_node_count));
    written += sizeof(_deleted_node_count);
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

void graph_t::load(std::istream& in) {
    //uint64_t written = 0;
    in.read((char*)&_max_node_id,sizeof(_max_node_id));
    in.read((char*)&_min_node_id,sizeof(_min_node_id));
    in.read((char*)&_node_count,sizeof(_node_count));
    in.read((char*)&_edge_count,sizeof(_edge_count));
    in.read((char*)&_path_count,sizeof(_path_count));
    in.read((char*)&_path_handle_next,sizeof(_path_handle_next));
    in.read((char*)&_deleted_node_count,sizeof(_deleted_node_count));
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
        uint64_t v;
        in.read((char*)&v,sizeof(uint64_t));
        path_name_map[string(k)] = v;
    }
}

}
