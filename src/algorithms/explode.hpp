#pragma once

#include <unordered_map>
#include <stack>
#include <map>
#include "odgi.hpp"
#include "handlegraph/path_position_handle_graph.hpp"
#include "weakly_connected_components.hpp"
#include "progress.hpp"

namespace odgi {
    namespace algorithms {

        /* This implementation has been taken from by: https://github.com/vgteam/vg */

        void add_connecting_edges_to_subgraph(const graph_t &source, graph_t &subgraph);

        void expand_subgraph_by_steps(const graph_t &source, graph_t &subgraph, const uint64_t &steps,
                                      bool forward_only);

        void add_full_paths_to_component(const graph_t &source, graph_t &component);
    }
}
