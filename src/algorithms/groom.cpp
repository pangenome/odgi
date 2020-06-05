/**
 * \file groom.cpp
 *
 * Defines an algorithm to remove spurious inverting links from the graph
 * by exploring the graph from the orientation supported by the most paths.
 */

#include "groom.hpp"

namespace odgi {
namespace algorithms {


void groom(handlegraph::MutablePathDeletableHandleGraph& source,
           handlegraph::MutablePathDeletableHandleGraph& target) {
           //bool use_heads, bool use_tails, bool show_progress) {

    //std::vector<handle_t> sorted = topological_order(&source, use_heads, use_tails, show_progress);
    /*
    ska::flat_hash_map<handlegraph::nid_t, std::pair<handlegraph::nid_t, bool>> splits
        = split_strands(&source, &target);
    */
    // the graph must have a compacted ID space
    //source.optimize(); if needed?
    /*
    dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,256,16> > seen, flipped;
    for (uint64_t i = 0; i < source.get_node_count(); ++i) {
        assert(graph.has_node(i+1));
        seen.push_back(0);
        flipped.push_back(0);
    }
    while (source.get_node_count() > target.get_node_count()) {
        uint64_t seed_rank = seen.select0(0);
        seen[seed_rank] = 1;
        handle_t seed = graph.get_handle(seed_rank+1);
        // start traversing from here
        
    }
    */

    // This (s) is our set of oriented nodes.
    //dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,256,16> > s;
    uint64_t min_handle_rank = 0;
    uint64_t max_handle_rank = 0;
    source.for_each_handle(
        [&](const handle_t& found) {
            uint64_t handle_rank = number_bool_packing::unpack_number(found);
            min_handle_rank = std::min(min_handle_rank, handle_rank);
            max_handle_rank = std::max(max_handle_rank, handle_rank);
        });

    // Start with the heads of the graph.
    // We could also just use the first node of the graph.
    std::vector<handle_t> seeds;
    bool use_heads = true;
    bool use_tails = false;
    if (use_heads) {
        seeds = head_nodes(&source);
    } else if (use_tails) {
        seeds = tail_nodes(&source);
    } else {
        handle_t min_handle = number_bool_packing::pack(min_handle_rank, false);
        seeds = { min_handle };
    }

    // We need to keep track of the nodes we haven't visited to seed subsequent
    // runs of the BFS
    dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,256,16> > unvisited, flipped;
    for (uint64_t i = 0; i <= max_handle_rank; ++i) {
        unvisited.push_back(1);
        flipped.push_back(0);
    }

    uint64_t prev_max_root = 0;
    uint64_t prev_max_length = 0;
    
    while (unvisited.rank1(unvisited.size())!=0) {

        bfs(source,
            [&source,&target,&unvisited,&flipped]
            (const handle_t& h, const uint64_t& r, const uint64_t& l, const uint64_t& d) {
                target.create_handle(source.get_sequence(h), source.get_id(h));
                uint64_t i = number_bool_packing::unpack_number(h);
                unvisited.set(i, 0);
                flipped.set(i, source.get_is_reverse(h));
            },
            [&unvisited](const handle_t& h) {
                uint64_t i = number_bool_packing::unpack_number(h);
                return unvisited.at(i)==0;
            },
            [](const handle_t& l, const handle_t& h) {
                // add the edge to the graph if it's not there yet
                //uint64_t i = number_bool_packing::unpack_number(h);
                //unvisited.set(i, 0);
                return false;
            },
            [](void) { return false; },
            seeds,
            { },
            false); // don't use bidirectional search
        // get another seed
        if (unvisited.rank1(unvisited.size())!=0) {
            uint64_t i = unvisited.select1(0);
            handle_t h = number_bool_packing::pack(i, false);
            seeds = { h };
        }
    }

    // add the edges
    source.for_each_edge(
        [&](const edge_t& edge) {
            //source.get_id(edge.first),
            bool from_flipped = flipped[number_bool_packing::unpack_number(edge.first)];
            handle_t from = target.get_handle(
                source.get_id(edge.first),
                source.get_is_reverse(edge.first)^from_flipped);
            bool to_flipped = flipped[number_bool_packing::unpack_number(edge.second)];
            handle_t to = target.get_handle(
                source.get_id(edge.second),
                source.get_is_reverse(edge.second)^to_flipped);
            target.create_edge(from, to);
        });
    // now add the paths back in
    source.for_each_path_handle(
        [&](const path_handle_t& path) {
            auto into = target.create_path_handle(source.get_path_name(path));
            source.for_each_step_in_path(
                path,
                [&](const step_handle_t& step) {
                    handle_t h = source.get_handle_of_step(step);
                    bool h_flipped = flipped[number_bool_packing::unpack_number(h)];
                    handle_t handle = target.get_handle(source.get_id(h),
                                                        source.get_is_reverse(h)^h_flipped);
                    target.append_step(into, handle);
                });
        });
    
    //std::cerr << "order size " << order.size() << " graph size " << g.get_node_count() << std::endl;
    //assert(order.size() == g.get_node_count());
    /*
    std::sort(order_raw.begin(), order_raw.end(),
              [](const bfs_state_t& a,
                 const bfs_state_t& b) {
                  return a.root < b.root || a.root == b.root && a.length < b.length;
              });

    std::vector<handle_t> order;
    for (auto& o : order_raw) order.push_back(o.handle);
    */
    //return order;

}


}
}

