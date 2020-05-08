/**
 * \file unchop.cpp
 *
 * Defines an algorithm to join adjacent handles.
 */

#include "unchop.hpp"

#include <unordered_set>
#include <list>
#include <set>
#include <iostream>
#include <sstream>

namespace odgi {
namespace algorithms {

/// Concatenates the nodes into a new node with the same external linkage as
/// the provided component. All handles must be in left to right order and a
/// consistent orientation. All paths present must run all the way through the
/// run of nodes from start to end or end to start.
///
/// Returns the handle to the newly created node.
///
/// After calling this on a vg::VG, paths will be invalid until
/// Paths::compact_ranks() is called.
handle_t concat_nodes(handlegraph::MutablePathDeletableHandleGraph& graph, const std::vector<handle_t>& nodes) {

    // Make sure we have at least 2 nodes
    assert(!nodes.empty() && nodes.front() != nodes.back());

    /*
#ifdef debug
    std::cerr << "Paths before concat: " << std::endl;
   
    graph.for_each_path_handle([&](const path_handle_t p) {
        std::cerr << graph.get_path_name(p) << ": ";
        for (auto h : graph.scan_path(p)) {
            std::cerr << graph.get_id(h) << (graph.get_is_reverse(h) ? '-' : '+') << " ";
        }
        std::cerr << std::endl;
    });
   
#endif
    */

    // We also require no edges enter or leave the run of nodes, but we can't check that now.
    
    // Make the new node
    handle_t new_node;
    {
        std::stringstream ss;
        for (auto& n : nodes) {
            ss << graph.get_sequence(n);
        }
        
        new_node = graph.create_handle(ss.str());
    }
    
#ifdef debug
    std::cerr << "Concatenating ";
    for (auto& n : nodes) {
        std::cerr << graph.get_id(n) << (graph.get_is_reverse(n) ? "-" : "+") << " ";
    }
    std::cerr << "into " << graph.get_id(new_node) << "+" << std::endl;
#endif
    
    // We should be able to rely on our handle graph to deduplicate edges, but see https://github.com/vgteam/libbdsg/issues/39
    // So we deduplicate ourselves.
    
    // Find all the neighbors. Make sure to translate edges to the other end of
    // the run, or self loops.
    std::unordered_set<handle_t> left_neighbors;
    graph.follow_edges(nodes.front(), true, [&](const handle_t& left_neighbor) {
        if (left_neighbor == nodes.back()) {
            // Loop back to the end
            left_neighbors.insert(new_node);
        } else if (left_neighbor == graph.flip(nodes.front())) {
            // Loop back to the front
            left_neighbors.insert(graph.flip(new_node));
        } else {
            // Normal edge
            left_neighbors.insert(left_neighbor);
        }
    });
    
    std::unordered_set<handle_t> right_neighbors;
    graph.follow_edges(nodes.back(), false, [&](const handle_t& right_neighbor) {
        if (right_neighbor == nodes.front()) {
            // Loop back to the front.
            // We will have seen it from the other side, so ignore it here.
        } else if (right_neighbor == graph.flip(nodes.back())) {
            // Loop back to the end
            right_neighbors.insert(graph.flip(new_node));
        } else {
            // Normal edge
            right_neighbors.insert(right_neighbor);
        }
    });
    
    // Make all the edges, now that we can't interfere with edge listing
    for (auto& n : left_neighbors) {
#ifdef debug
        std::cerr << "Creating edge " << graph.get_id(n) << (graph.get_is_reverse(n) ? "-" : "+") << " -> "
            <<  graph.get_id(new_node) << (graph.get_is_reverse(new_node) ? "-" : "+") << std::endl;
#endif
        graph.create_edge(n, new_node);
    }
    for (auto& n : right_neighbors) {
    
#ifdef debug
        std::cerr << "Creating edge " << graph.get_id(new_node) << (graph.get_is_reverse(new_node) ? "-" : "+") << " -> "
            <<  graph.get_id(n) << (graph.get_is_reverse(n) ? "-" : "+") << std::endl;
#endif
    
        graph.create_edge(new_node, n);
    }
    
    {
        // Collect the first and last visits along paths. TODO: this requires
        // walking each path all the way through.
        //
        // This contains the first and last handles in path orientation, and a flag
        // for if the path runs along the reverse strand of our run of nodes.
        std::vector<std::tuple<step_handle_t, step_handle_t, bool>> ranges_to_rewrite;
        
        graph.for_each_step_on_handle(nodes.front(), [&](const step_handle_t& front_step) {
            auto path = graph.get_path_handle_of_step(front_step);
#ifdef debug
            std::cerr << "Consider path " << graph.get_path_name(path) << std::endl;
#endif
        
            // If we don't get the same oriented node as where this step is
            // stepping, we must be stepping on the other orientation.
            bool runs_reverse = (graph.get_handle_of_step(front_step) != nodes.front());
       
            step_handle_t back_step = front_step;
            
            while(graph.get_handle_of_step(back_step) != (runs_reverse ? graph.flip(nodes.back()) : nodes.back())) {
                // Until we find the step on the path that visits the last node in our run.
                // Go along the path towards where our last node should be, in our forward orientation.
                back_step = runs_reverse ? graph.get_previous_step(back_step) : graph.get_next_step(back_step);
            }
            
            // Now we can record the range to rewrite
            // Make sure to put it into path-forward order
            if (runs_reverse) {
#ifdef debug
                std::cerr << "\tGoing to rewrite between " << graph.get_id(graph.get_handle_of_step(front_step)) << " and " << graph.get_id(graph.get_handle_of_step(back_step)) << " backward" << std::endl;
#endif
                ranges_to_rewrite.emplace_back(back_step, front_step, true);
            } else {
            
#ifdef debug
                std::cerr << "\tGoing to rewrite between " << graph.get_id(graph.get_handle_of_step(front_step)) << " and " << graph.get_id(graph.get_handle_of_step(back_step)) << std::endl;
#endif
                ranges_to_rewrite.emplace_back(front_step, back_step, false);
            }
        });

        uint64_t i = 0;
        for (auto& range : ranges_to_rewrite) {
            // Rewrite each range to visit the new node in the appropriate orientation instead of whatever it did before
            // Make sure to advance the end of the range because rewrite is end-exclusive (to allow insert).
            graph.rewrite_segment(std::get<0>(range), std::get<1>(range), {std::get<2>(range) ? graph.flip(new_node) : new_node});
        }
    }

    // Destroy all the old edges
    // We know they only exist to the left and right neighbors, and along the run
    for (auto& n : left_neighbors) {
        graph.destroy_edge(n, nodes.front());
    }
    for (auto& n : right_neighbors) {
        graph.destroy_edge(nodes.back(), n);
    }
    auto it = nodes.begin();
    auto next_it = it;
    ++next_it;
    while (next_it != nodes.end()) {
        graph.destroy_edge(*it, *next_it);
        it = next_it;
        ++next_it;
    }
    
    for (auto& n : nodes) {
        // Destroy all the old nodes
#ifdef debug
        std::cerr << "Destroying node " << graph.get_id(n) << std::endl;
#endif
        graph.destroy_handle(n);
    }

    /*
#ifdef debug
    std::cerr << "Paths after concat: " << std::endl;
   
    graph.for_each_path_handle([&](const path_handle_t p) {
        std::cerr << graph.get_path_name(p) << ": ";
        for (auto h : graph.scan_path(p)) {
            std::cerr << graph.get_id(h) << (graph.get_is_reverse(h) ? '-' : '+') << " ";
        }
        std::cerr << std::endl;
    });
   
#endif
    */

    // Return the new handle we merged to.
    return new_node;
}

void unchop(handlegraph::MutablePathDeletableHandleGraph& graph) {
#ifdef debug
    std::cerr << "Running unchop" << std::endl;
#endif
    for (auto& comp : simple_components(graph, 2)) {
#ifdef debug
        std::cerr << "Unchop " << comp.size() << " nodes together" << std::endl;
#endif
        concat_nodes(graph, comp);
    }
}


}
}

