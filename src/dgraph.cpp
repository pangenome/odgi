//
//  dgraph.cpp
//  
//
//  Created by Jordan Eizenga on 12/11/18.
//

#include "dgraph.hpp"

namespace dankgraph {
    
    const double SuccinctDynamicSequenceGraph::defrag_factor = .2;
    
    SuccinctDynamicSequenceGraph::SuccinctDynamicSequenceGraph() {
        
    }
    
    SuccinctDynamicSequenceGraph::~SuccinctDynamicSequenceGraph() {
        
    }
    
    handle_t SuccinctDynamicSequenceGraph::create_handle(const std::string& sequence) {
        return create_handle(sequence, max_id + 1);
    }
    
    handle_t SuccinctDynamicSequenceGraph::create_handle(const std::string& sequence, const id_t& id) {
        
        // TODO: don't duplicate an existing node
        
        // no edges yet, null pointer for linked list
        graph_iv.append(0);
        graph_iv.append(0);
        // record the sequence interval
        graph_iv.append(seq_iv.size());
        graph_iv.append(sequence.size());
        
        // encode the sequence interval
        for (size_t i = 0; i < sequence.size(); i++) {
            seq_iv.append(encode_nucleotide(sequence[i]));
            boundary_bv.append(i ? 0 : 1);
        }
        
        // expand the ID vector's dimensions so it can handle the full ID interval
        if (id_to_graph_iv.empty()) {
            id_to_graph_iv.append_back(0);
        }
        else {
            for (int64_t i = id; i < min_id; i++) {
                id_to_graph_iv.append_front(0);
            }
            for (int64_t i = id_to_graph_iv.size(); i <= id - min_id; i++) {
                id_to_graph_iv.append_back(0);
            }
        }
        
        // update the min and max ID
        max_id = std::max(id, max_id);
        min_id = std::max(id, min_id);
        
        // record the mapping of the ID to the graph record
        id_to_graph_iv.set(id - min_id, graph_iv.size() / GRAPH_RECORD_SIZE);
        
        return get_handle(id);
    }
    
    void SuccinctDynamicSequenceGraph::create_edge(const handle_t& left, const handle_t& right) {
        
        // look for the edge
        bool add_edge = follow_edges(left, false, [&](const handle_t& next) {
            return next == right;
        });
        
        // don't duplicate it
        if (!add_edge) {
            return;
        }
        
        // get the location of the edge list pointer in the graph vector
        size_t g_iv_left = graph_iv_index(left) + (get_is_reverse(left) ?
                                                   GRAPH_START_EDGES_OFFSET :
                                                   GRAPH_END_EDGES_OFFSET);
        size_t g_iv_right = graph_iv_index(left) + (get_is_reverse(left) ?
                                                    GRAPH_END_EDGES_OFFSET :
                                                    GRAPH_START_EDGES_OFFSET);
        
        // add a new linked list node pointing to the rest of the list
        edge_lists_iv.append(encode_edge_target(right));
        edge_lists_iv.append(graph_iv.get(g_iv_left));
        // make this new node the head
        graph_iv.set(g_iv_left, edge_lists_iv.size() / EDGE_RECORD_SIZE);
        
        // don't double add a reversing self edge
        if (g_iv_left == g_iv_right) {
            return;
        }
        
        // add a new linked list node pointing to the rest of the list
        edge_lists_iv.append(encode_edge_target(flip(left)));
        edge_lists_iv.append(graph_iv.get(g_iv_right));
        // make this new node the head
        graph_iv.set(g_iv_right, edge_lists_iv.size() / EDGE_RECORD_SIZE);
    }
    
    handle_t SuccinctDynamicSequenceGraph::get_handle(const id_t& node_id, bool is_reverse) const {
        return EasyHandlePacking::pack(node_id, is_reverse);
    }
    
    id_t SuccinctDynamicSequenceGraph::get_id(const handle_t& handle) const {
        return EasyHandlePacking::unpack_number(handle);
    }
    
    bool SuccinctDynamicSequenceGraph::get_is_reverse(const handle_t& handle) const {
        return EasyHandlePacking::unpack_bit(handle);;
    }
    
    handle_t SuccinctDynamicSequenceGraph::flip(const handle_t& handle) const {
        return EasyHandlePacking::toggle_bit(handle);
    }
    
    size_t SuccinctDynamicSequenceGraph::get_length(const handle_t& handle) const {
        return graph_iv.get(graph_iv_index(handle) + GRAPH_SEQ_LENGTH_OFFSET);
    }
    
    std::string SuccinctDynamicSequenceGraph::get_sequence(const handle_t& handle) const {
        size_t g_iv_index = graph_iv_index(handle);
        size_t seq_start = graph_iv.get(g_iv_index + GRAPH_SEQ_START_OFFSET);
        size_t seq_len = graph_iv.get(g_iv_index + GRAPH_SEQ_LENGTH_OFFSET);
        std::string seq(seq_len, 'A');
        for (size_t i = 0; i < seq_len; i++) {
            seq[i] = decode_nucleotide(seq_iv.get(seq_start + i));
        }
        return get_is_reverse(handle) ? reverse_complement(seq) : seq;
    }
    
    void SuccinctDynamicSequenceGraph::swap_handles(const handle_t& a, const handle_t& b) {
        size_t g_iv_index_a = graph_iv_index(a);
        size_t g_iv_index_b = graph_iv_index(b);
        
        int64_t id_a = get_id(a);
        int64_t id_b = get_id(b);
        
        for (size_t i = 0; i < GRAPH_RECORD_SIZE; i++) {
            uint64_t val = graph_iv.get(g_iv_index_a + i);
            graph_iv.set(g_iv_index_a + i, graph_iv.get(g_iv_index_b + i));
            graph_iv.set(g_iv_index_b + i, val);
        }
        
        id_to_graph_iv.set(id_a, g_iv_index_b);
        id_to_graph_iv.set(id_b, g_iv_index_a);
    }
    
    bool SuccinctDynamicSequenceGraph::follow_edges(const handle_t& handle, bool go_left,
                                                    const std::function<bool(const handle_t&)>& iteratee) const {
        // toward start = true, toward end = false
        bool direction = get_is_reverse(handle) != go_left;
        // get the head of the linked list from the graph vector
        size_t edge_idx = graph_iv.get(graph_iv_index(handle)
                                       + (direction ? GRAPH_START_EDGES_OFFSET : GRAPH_END_EDGES_OFFSET));
        // traverse the linked list as long as directed
        bool keep_going = true;
        while (edge_idx && keep_going) {
            
            handle_t edge_target = decode_edge_target(get_edge_target(edge_idx));
            if (go_left) {
                // match the orientation encoding
                edge_target = flip(edge_target);
            }
            
            keep_going = iteratee(edge_target);
            edge_idx = get_next_edge_index(edge_idx);
        }
        
        return keep_going;
    }
    
    size_t SuccinctDynamicSequenceGraph::node_size(void) const {
        return graph_iv.size() / GRAPH_RECORD_SIZE;
    }
    
    id_t SuccinctDynamicSequenceGraph::min_node_id(void) const {
        return min_id;
    }
    
    id_t SuccinctDynamicSequenceGraph::max_node_id(void) const {
        return max_id;
    }
    
    void SuccinctDynamicSequenceGraph::for_each_handle(const std::function<bool(const handle_t&)>& iteratee,
                                                       bool parallel) const {
        
        size_t num_id_positions = id_to_graph_iv.size();
        if (parallel) {
            // TODO: add OMP pragma back in
            // TODO: would task based parallelism be better?
//#pragma omp parallel for
            for (size_t i = 0; i < num_id_positions; i++) {
                if (id_to_graph_iv.get(i)) {
                    iteratee(get_handle(i + min_id));
                }
            }
        }
        else {
            for (size_t i = 0; i < num_id_positions; i++) {
                if (id_to_graph_iv.get(i)) {
                    iteratee(get_handle(i + min_id));
                }
            }
        }
        
    }
    
    handle_t SuccinctDynamicSequenceGraph::apply_orientation(const handle_t& handle) {
        
        if (get_is_reverse(handle)) {
            size_t g_iv_idx = graph_iv_index(handle);
            size_t seq_start = graph_iv.get(g_iv_idx + GRAPH_SEQ_START_OFFSET);
            size_t seq_len = graph_iv.get(g_iv_idx + GRAPH_SEQ_LENGTH_OFFSET);
            
            // reverse complement the sequence in place
            for (size_t i = 0; i < seq_len / 2; i++) {
                size_t j = seq_start + seq_len - i - 1;
                size_t k = seq_start + i;
                uint64_t base = seq_iv.get(k);
                seq_iv.set(k, complement_encoded_nucleotide(seq_iv.get(j)));
                seq_iv.set(j, complement_encoded_nucleotide(base));
            }
            if (seq_len % 2) {
                size_t j = seq_start + seq_len / 2;
                seq_iv.set(j, complement_encoded_nucleotide(seq_iv.get(j)));
            }
            
            // the ID is preserved, we just need to need to return a forward version
            return flip(handle);
        }
        else {
            // it's already the way we want it
            return handle;
        }
    }
    
    std::vector<handle_t> SuccinctDynamicSequenceGraph::divide_handle(const handle_t& handle,
                                                                      const std::vector<size_t>& offsets) {
        
        // put the offsets in forward orientation to simplify subsequent steps
        std::vector<size_t> forward_offsets = offsets;
        size_t node_length = get_length(handle);
        if (get_is_reverse(handle)) {
            for (size_t& off : forward_offsets) {
                off = node_length - off;
            }
        }
        
        // we will also build the return value in forward orientation
        handle_t forward_handle = get_is_reverse(handle) ? flip(handle) : handle;
        std::vector<handle_t> return_val{forward_handle};
        size_t g_iv_idx = graph_iv_index(forward_handle);
        
        // offsets in the sequence vector will be measured relative to the first position of
        // the current handle
        size_t first_start = graph_iv.get(g_iv_idx + GRAPH_SEQ_START_OFFSET);
        
        // we record the the edges out of this node so they can be transferred onto the final
        // node in the split
        size_t end_edges_idx = graph_iv.get(g_iv_idx + GRAPH_END_EDGES_OFFSET);
        
        // init trackers for the previous iteration
        size_t last_offset = 0;
        id_t prev_id = get_id(forward_handle);
        for (const size_t& off : offsets) {
            // make a new node record:
            // ID
            graph_iv.append(max_id + 1);
            // start edges
            graph_iv.append(0);
            // end edges
            graph_iv.append(0);
            // seq start
            graph_iv.append(first_start + off);
            boundary_bv.set(first_start + off, 1);
            // set len (just a placeholder for now, will be set in next iteration)
            graph_iv.append(0);
            // record the mapping of the ID to the graph record in the final position
            id_to_graph_iv.append_back(graph_iv.size() / GRAPH_RECORD_SIZE);
            // add the new handle to the return vector
            return_val.push_back(get_handle(max_id + 1, false));
            
            // now let's do what we still need to on the previous node
            // set the previous node's length based on the current offset
            graph_iv.set(g_iv_idx + GRAPH_SEQ_LENGTH_OFFSET, off - last_offset);
            // create an edge to the new node
            edge_lists_iv.append(encode_edge_target(get_handle(max_id + 1)));
            edge_lists_iv.append(0);
            // add the edge onto the previous node
            graph_iv.set(g_iv_idx + GRAPH_END_EDGES_OFFSET, edge_lists_iv.size() / EDGE_RECORD_SIZE);
            
            // move the pointer to the node we just made
            g_iv_idx = graph_iv.size() / GRAPH_RECORD_SIZE;
            
            // create an edge backward
            edge_lists_iv.append(encode_edge_target(get_handle(prev_id, true)));
            edge_lists_iv.append(0);
            // add the edge backwards to the current node
            graph_iv.set(g_iv_idx + GRAPH_START_EDGES_OFFSET, edge_lists_iv.size() / EDGE_RECORD_SIZE);
            
            // update the offset and ID trackers
            prev_id = max_id + 1;
            last_offset = off;
            max_id++;
        }
        // set final node's length to the remaining sequence
        graph_iv.set(g_iv_idx + GRAPH_SEQ_LENGTH_OFFSET, node_length - last_offset);
        
        // point the final node's end edges to the original node's end edges
        graph_iv.set(g_iv_idx + GRAPH_END_EDGES_OFFSET, end_edges_idx);
        
        if (get_is_reverse(handle)) {
            // reverse the vector to the orientation of the input handle
            std::reverse(return_val.begin(), return_val.end());
            for (handle_t& ret_handle : return_val) {
                ret_handle = flip(ret_handle);
            }
        }
        
        return return_val;
    }
    
    void SuccinctDynamicSequenceGraph::destroy_handle(const handle_t& handle) {
        
        // remove the reference to the node
        id_to_graph_iv.set(get_id(handle), 0);
        
        // remove the back-references to the edges
        follow_edges(handle, false, [&](const handle_t& next) {
            remove_edge_reference(flip(next), flip(handle));
            // we don't actually bother removing the reference, but we will also consider
            // the edge on the deleting node to be deleted
            deleted_edge_records++;
            return true;
        });
        follow_edges(handle, true, [&](const handle_t& prev) {
            remove_edge_reference(prev, handle);
            // we don't actually bother removing the reference, but we will also consider
            // the edge on the deleting node to be deleted
            deleted_edge_records++;
            return true;
        });
        
        deleted_node_records++;
        
        // maybe reallocate to address fragmentation
        defragment();
    }
    
    void SuccinctDynamicSequenceGraph::remove_edge_reference(const handle_t& on, const handle_t& to) {
        
        handle_t looking_for = get_is_reverse(on) ? flip(to) : to;
        
        size_t g_iv_idx = graph_iv_index(on) + (get_is_reverse(on)
                                                ? GRAPH_END_EDGES_OFFSET
                                                : GRAPH_START_EDGES_OFFSET);
        size_t edge_list_idx = graph_iv.get(g_iv_idx);
        
        if (decode_edge_target(get_edge_target(edge_list_idx)) == looking_for) {
            // the edge back to the deleting node is the first in the list, so we need
            // to update the head
            graph_iv.set(g_iv_idx, get_next_edge_index(edge_list_idx));
        }
        else {
            // we need to traverse down the list and to find the edge back
            size_t prev_edge_list_idx = edge_list_idx;
            edge_list_idx = get_next_edge_index(edge_list_idx);
            while (decode_edge_target(get_edge_target(edge_list_idx)) != looking_for) {
                prev_edge_list_idx = edge_list_idx;
                edge_list_idx = get_next_edge_index(edge_list_idx);
            }
            // skip over this edge in this linked list
            edge_lists_iv.set((prev_edge_list_idx - 1) * EDGE_RECORD_SIZE + EDGE_NEXT_OFFSET,
                              get_next_edge_index(edge_list_idx));
        }
        deleted_edge_records++;
    }
    
    void SuccinctDynamicSequenceGraph::destroy_edge(const handle_t& left, const handle_t& right) {
        remove_edge_reference(left, right);
        remove_edge_reference(flip(right), flip(left));
        defragment();
    }
    
    void SuccinctDynamicSequenceGraph::defragment(void) {
        if (deleted_node_records > defrag_factor * (graph_iv.size() / GRAPH_RECORD_SIZE)) {
            // adjust the start
            while (id_to_graph_iv.empty() ? false : id_to_graph_iv.get(0) == 0) {
                id_to_graph_iv.pop_front();
                min_id++;
            }
            // adjust the end
            while (id_to_graph_iv.empty() ? false : id_to_graph_iv.get(0) == 0) {
                id_to_graph_iv.pop_back();
            }
            
            SuccinctDynamicVector new_graph_iv;
            
            for (size_t i = 0; i < id_to_graph_iv.size(); i++) {
                size_t g_iv_idx = id_to_graph_iv.get(i);
                if (g_iv_idx) {
                    // this node still exists, create a new copy
                    new_graph_iv.append(graph_iv.get(g_iv_idx + GRAPH_START_EDGES_OFFSET));
                    new_graph_iv.append(graph_iv.get(g_iv_idx + GRAPH_END_EDGES_OFFSET));
                    new_graph_iv.append(graph_iv.get(g_iv_idx + GRAPH_SEQ_START_OFFSET));
                    new_graph_iv.append(graph_iv.get(g_iv_idx + GRAPH_SEQ_LENGTH_OFFSET));
                    // update the pointer into graph_iv
                    id_to_graph_iv.set(i, new_graph_iv.size() / GRAPH_RECORD_SIZE);
                }
            }
            
            // replace graph_iv with the defragged copy
            graph_iv = new_graph_iv;
            deleted_node_records = 0;
        }
        
        // TODO: defrag the seq_iv?
        
        if (deleted_edge_records > defrag_factor * (edge_lists_iv.size() / EDGE_RECORD_SIZE)) {
            
            SuccinctDynamicVector new_edge_lists_iv;
            
            for (size_t i = 0; i < id_to_graph_iv.size(); i++) {
                size_t g_iv_idx = id_to_graph_iv.get(i);
                if (g_iv_idx) {
                    for (bool on_start : {true, false}) {
                        // are we doing the start or the end edges?
                        size_t which_edge_offset = on_start ? GRAPH_START_EDGES_OFFSET : GRAPH_END_EDGES_OFFSET;
                        
                        size_t edge_list_idx = graph_iv.get(g_iv_idx + which_edge_offset);
                        if (edge_list_idx) {
                            // add a new edge record
                            new_edge_lists_iv.append(get_edge_target(edge_list_idx));
                            new_edge_lists_iv.append(0);
                            // point the graph vector at this new edge list
                            graph_iv.set(g_iv_idx + which_edge_offset, new_edge_lists_iv.size() / EDGE_RECORD_SIZE);
                            
                            edge_list_idx = get_next_edge_index(edge_list_idx);
                            while (edge_list_idx) {
                                // add a new edge record
                                new_edge_lists_iv.append(get_edge_target(edge_list_idx));
                                new_edge_lists_iv.append(0);
                                // point the previous link at this one
                                new_edge_lists_iv.set(new_edge_lists_iv.size() - 2 * EDGE_RECORD_SIZE + EDGE_NEXT_OFFSET,
                                                      new_edge_lists_iv.size() / EDGE_RECORD_SIZE);
                                
                                edge_list_idx = get_next_edge_index(edge_list_idx);
                            }
                        }
                    }
                }
            }
            
            edge_lists_iv = new_edge_lists_iv;
            
            deleted_edge_records = 0;
        }
    }
    
    void SuccinctDynamicSequenceGraph::clear(void) {
        graph_iv.clear();
        edge_lists_iv.clear();
        id_to_graph_iv.clear();
        seq_iv.clear();
        boundary_bv.clear();
        min_id = std::numeric_limits<id_t>::max();
        max_id = 0;
        deleted_edge_records = 0;
        deleted_node_records = 0;
    }
}
