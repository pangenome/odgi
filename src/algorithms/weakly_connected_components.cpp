#include "weakly_connected_components.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

std::vector<ska::flat_hash_set<handlegraph::nid_t>> weakly_connected_components(const HandleGraph* graph) {
    std::vector<ska::flat_hash_set<handlegraph::nid_t>> to_return;
    
    // This only holds locally forward handles
    ska::flat_hash_set<handle_t> traversed;
    
    graph->for_each_handle([&](const handle_t& handle) {
        
        // Only think about it in the forward orientation
        auto forward = graph->forward(handle);
        
        if (traversed.count(forward)) {
            // Already have this node, so don't start a search from it.
            return;
        }
        
        // The stack only holds locally forward handles
        std::vector<handle_t> stack{forward};
        to_return.emplace_back();
        while (!stack.empty()) {
            handle_t here = stack.back();
            stack.pop_back();
            
            traversed.insert(here);
            to_return.back().insert(graph->get_id(here));
            
            // We have a function to handle all connected handles
            auto handle_other = [&](const handle_t& other) {
                // Again, make it forward
                auto other_forward = graph->forward(other);
                
                if (!traversed.count(other_forward)) {
                    stack.push_back(other_forward);
                }
            };
            
            // Look at edges in both directions
            graph->follow_edges(here, false, handle_other);
            graph->follow_edges(here, true, handle_other);
            
        }
    });
    return to_return;
}

std::vector<std::vector<handlegraph::handle_t>> weakly_connected_component_vectors(const HandleGraph* graph) {
    std::vector<std::vector<handlegraph::handle_t>> components;
    for (auto& component : weakly_connected_components(graph)) {
        components.emplace_back();
        auto& v = components.back();
        for (auto& id : component) {
            v.push_back(graph->get_handle(id));
        }
        std::sort(v.begin(), v.end(),
                  [](const handle_t& a,
                     const handle_t& b) {
                      return as_integer(a) < as_integer(b);
                  });
    }
    return components;
}

std::vector<std::pair<ska::flat_hash_set<handlegraph::nid_t>, std::vector<handle_t>>> weakly_connected_components_with_tips(const HandleGraph* graph) {
    // TODO: deduplicate with above
    
    std::vector<std::pair<ska::flat_hash_set<handlegraph::nid_t>, std::vector<handle_t>>> to_return;
    
    // This only holds locally forward handles
    ska::flat_hash_set<handle_t> traversed;
    
    graph->for_each_handle([&](const handle_t& handle) {
        
        // Only think about it in the forward orientation
        auto forward = graph->forward(handle);
        
        if (traversed.count(forward)) {
            // Already have this node, so don't start a search from it.
            return;
        }
        
        // The stack only holds locally forward handles
        std::vector<handle_t> stack{forward};
        to_return.emplace_back();
        while (!stack.empty()) {
            handle_t here = stack.back();
            stack.pop_back();
            
            traversed.insert(here);
            to_return.back().first.insert(graph->get_id(here));
            
            // We have a counter for the number of edges we saw.
            // If it is 0 after traversing edges we know we have a tip.
            size_t edge_counter = 0;
            
            // We have a function to handle all connected handles
            auto handle_other = [&](const handle_t& other) {
                // Again, make it forward
                auto other_forward = graph->forward(other);
                
                if (!traversed.count(other_forward)) {
                    stack.push_back(other_forward);
                }
                
                edge_counter++;
            };
            
            // Look at edges in both directions
            graph->follow_edges(here, false, handle_other);
            if (edge_counter == 0) {
                // This is a tail node. Put it in reverse as a tip.
                to_return.back().second.push_back(graph->flip(here));
            }
            
            edge_counter = 0;
            graph->follow_edges(here, true, handle_other);
            if (edge_counter == 0) {
                // This is a head node. Put it as a tip.
                to_return.back().second.push_back(here);
            }
            
        }
    });
    return to_return;
}

}
}
