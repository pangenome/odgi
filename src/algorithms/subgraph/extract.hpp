#ifndef ODGI_EXTRACT_H
#define ODGI_EXTRACT_H

#include <unordered_map>
#include <stack>
#include <map>
#include "odgi.hpp"
#include "handlegraph/path_position_handle_graph.hpp"
#include "weakly_connected_components.hpp"
#include "progress.hpp"

#include "src/algorithms/subgraph/region.hpp"

namespace odgi {
    namespace algorithms {
        /* Several functions were inspired by https://github.com/vgteam/vg */

        void add_full_paths_to_component(const graph_t &source, graph_t &component);

        /// add subpaths to the subgraph, providing a concatenation of subpaths that are discontiguous over the subgraph
        /// based on their order in the path position index provided by the source graph
        /// will clear any path found in both graphs before writing the new steps into it
        /// a suffix will be added to each path in the subgraph denoting its offset
        /// in the source graph (unless the subpath was not cut up at all)
        void add_subpaths_to_subgraph(const graph_t &source, graph_t &subgraph, uint64_t num_threads);

        void extract_path_range(const graph_t &source, path_handle_t path_handle, int64_t start, int64_t end,
                                graph_t &subgraph);

        void add_connecting_edges_to_subgraph(const graph_t &source, graph_t &subgraph);

        void expand_subgraph_by_length(const graph_t &source, graph_t &subgraph, const uint64_t &length,
                                       bool forward_only);

        void expand_subgraph_by_steps(const graph_t &source, graph_t &subgraph, const uint64_t &steps,
                                      bool forward_only);

        void extract_id_range(const graph_t& source, const nid_t& id1, const nid_t& id2, graph_t& subgraph);

    }
}
#endif //ODGI_EXTRACT_H
