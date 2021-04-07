/**
 * \file unroll.cpp
 *
 * Linear unrolling of the graph
 */

#include "unroll.hpp"

namespace odgi {
namespace algorithms {

void unroll(
    const handlegraph::PathHandleGraph& input,
    handlegraph::MutablePathDeletableHandleGraph& output) {

    // map from id in input to a vector of previous ranks and node ids in the new graph
    std::vector<std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>> node_map(input.get_node_count()+1);
    // add the base nodes
    input.for_each_handle(
        [&](const handle_t& handle) {
            auto id = input.get_id(handle);
            node_map[id].push_back({0, 0, id});
            auto new_handle = output.create_handle(input.get_sequence(handle));
            assert(output.get_id(new_handle) == id);
        });

    auto get_unrolled_handle =
        [&](const handle_t& handle,
            const uint64_t& count,
            const uint64_t& prev_count,
            const uint64_t& prev_id,
            uint64_t& assigned_level) {
            auto id = input.get_id(handle);
            auto& new_nodes = node_map[id];
            uint64_t i = 0;
            for (auto& n : new_nodes) {
                if (i == count && std::get<1>(n) == prev_count
                    && (std::get<0>(n) == 0 || std::get<0>(n) == prev_id)) {
                    break;
                }
                ++i;
            }
            if (i == new_nodes.size()) {
                assigned_level = i;
                // create a new node in the output graph
                auto new_handle = output.create_handle(input.get_sequence(input.get_handle(id)));
                new_nodes.push_back({prev_id, prev_count, output.get_id(new_handle)});
                return input.get_is_reverse(handle) ? output.flip(new_handle) : new_handle;
            } else {
                auto new_handle = output.get_handle(std::get<2>(new_nodes[i]));
                return input.get_is_reverse(handle) ? output.flip(new_handle) : new_handle;
            }
        };
            
    input.for_each_path_handle(
        [&](const path_handle_t& path) {
            ska::flat_hash_map<uint64_t, uint64_t> seen_count; // how many times we've seen the node in this path
            auto unrolled_path = output.create_path_handle(input.get_path_name(path));
            uint64_t prev_count = 0;
            uint64_t prev_id = 0;
            auto handle_step =
                [&](const step_handle_t& step) {
                    auto handle = input.get_handle_of_step(step);
                    auto id = input.get_id(handle);
                    // is it the first time we're here?
                    auto& count = seen_count[id];
                    uint64_t assigned_level = 0;
                    output.append_step(unrolled_path,
                                       get_unrolled_handle(
                                           handle, count, prev_count, prev_id, assigned_level));
                    prev_count = (count == 0 ? 0 : assigned_level); //std::max(assigned_level, count));
                    prev_id = id;
                    ++count;
                };

            input.for_each_step_in_path(
                path, handle_step);
        });

    // embed all paths in the graph to ensure validity
    output.for_each_path_handle(
        [&](const path_handle_t& path) {
            handle_t last;
            step_handle_t begin_step = output.path_begin(path);
            output.for_each_step_in_path(
                path,
                [&](const step_handle_t &step) {
                    handle_t h = output.get_handle_of_step(step);
                    if (step != begin_step) {
                        output.create_edge(last, h);
                    }
                    last = h;
                });
        });

    // now link all the paths together in the unrolled graph
    
    
    
}


}
}
