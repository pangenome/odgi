/**
 * \file chop.cpp
 *
 * Defines an algorithm to set a maximum node length by dividing the graph nodes
 */

#include "chop.hpp"

#include <vector>
#include <iostream>

namespace odgi {
namespace algorithms {

void chop(handlegraph::MutablePathDeletableHandleGraph& graph,
          const uint64_t& max_node_length) {
    std::vector<handle_t> to_chop;
    graph.for_each_handle([&](const handle_t& handle) {
            if (graph.get_length(handle) > max_node_length) {
                to_chop.push_back(handle);
            }
        });
    for (auto& handle : to_chop) {
        // get divide points
        uint64_t length = graph.get_length(handle);
        std::vector<size_t> offsets;
        for (uint64_t i = max_node_length; i < length; i+=max_node_length) {
            offsets.push_back(i);
        }
        graph.divide_handle(handle, offsets);
    }
}

}
}

