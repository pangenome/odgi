/**
 * \file chop.cpp
 *
 * Defines an algorithm to set a maximum node length by dividing the graph nodes
 */

#include "chop.hpp"

#include <vector>
#include <iostream>
#include <deps/ips4o/ips4o.hpp>

namespace odgi {
    namespace algorithms {

        void chop(handlegraph::MutablePathDeletableHandleGraph &graph,
                  const uint64_t &max_node_length, const uint64_t &nthreads, const bool &show_info) {
            std::vector<std::tuple<uint64_t, uint64_t, handle_t>> originalRank_inChoppedNodeRank_handle;
            std::vector<std::pair<uint64_t, handle_t>> originalRank_handleToChop;
            uint64_t rank = 0;
            graph.for_each_handle([&](const handle_t &handle) {
                if (graph.get_length(handle) > max_node_length) {
                    originalRank_handleToChop.push_back(std::make_pair(rank, handle));
                } else {
                    originalRank_inChoppedNodeRank_handle.push_back(std::make_tuple(rank, 0, handle));
                }

                rank++;
            });

            if (show_info) {
                std::cerr << "[odgi::chop] " << originalRank_handleToChop.size() << " node(s) to chop." << std::endl;
            }

            for (auto rank_handle : originalRank_handleToChop) {
                // get divide points
                uint64_t length = graph.get_length(rank_handle.second);
                std::vector<size_t> offsets;
                for (uint64_t i = max_node_length; i < length; i += max_node_length) {
                    offsets.push_back(i);
                }

                rank = 0;
                for (auto chopped_handle : graph.divide_handle(rank_handle.second, offsets)) {
                    originalRank_inChoppedNodeRank_handle.push_back(
                            std::make_tuple(rank_handle.first, rank, chopped_handle));

                    rank++;
                }
            }

            ips4o::parallel::sort(
                    originalRank_inChoppedNodeRank_handle.begin(), originalRank_inChoppedNodeRank_handle.end(),
                    std::less<>(), nthreads
            );

            std::vector<handle_t> new_handles;
            for (auto x_y_z : originalRank_inChoppedNodeRank_handle) {
                new_handles.push_back(std::get<2>(x_y_z));
            }

            graph.apply_ordering(new_handles, true);
        }

    }
}

