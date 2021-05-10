#pragma once

#include <unordered_map>
#include <stack>
#include <map>
#include <sdsl/bit_vectors.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
#include <handlegraph/util.hpp>
#include "weakly_connected_components.hpp"
#include <deps/ips4o/ips4o.hpp>
#include "XoshiroCpp.hpp"
#include "progress.hpp"

namespace odgi {
namespace algorithms {

        using namespace handlegraph;

        /* This implementation has been inspired by: https://github.com/jltsiren/gbwtgraph */

        constexpr size_t PATH_COVER_DEFAULT_N = 16;
        constexpr size_t PATH_COVER_DEFAULT_K = 2;

        /*
          Determine whether the given component is acyclic in a nice way. In such graphs,
          when we start from nodes with indegree 0 in forward orientation, we reach each node
          in a single orientation and find no cycles. Return the head nodes when the component
          passes the tests or an empty vector otherwise.
          Ignores node ids that are not present in the graph.
        */
        ska::flat_hash_set<handlegraph::nid_t>
        is_nice_and_acyclic(const HandleGraph &graph, const ska::flat_hash_set<handlegraph::nid_t> &component);

        /*
          Find a path cover of the graph with num_paths_per_component paths per component, adding the generated paths in
          the graph. The path cover is built greedily. Each time we extend a path, we choose the extension,
          where the coverage of the node_window_size >= 2 node window is the lowest. Note that this is a maximum
          coverage algorithm that tries to maximize the number of windows covered by a fixed number
          of paths.

          This algorithm has been inspired by:
            Ghaffaari and Marschall: Fully-sensitive seed finding in sequence graphs using a
            hybrid index. Bioinformatics, 2019.

          Because the graph here may have cycles and orientation flips, some changes were necessary:
            - If the component is not a DAG, we start from an arbitrary node with minimal
              coverage and extend in both directions.
            - In a DAG, we start from the head node with the lowest coverage so far.
            - We stop when the length of the path reaches the size of the component.
            - The length of the window is in nodes instead of base pairs. We expect a sparse graph,
              where the nodes between variants are long.
            - If the component is not a DAG and the path is shorter than node_window_size - 1 nodes, we consider
              the coverage of individual nodes instead of windows.
            - When determining window coverage, we consider the window equivalent to its reverse
              complement.
        */
        void path_cover(handlegraph::MutablePathDeletableHandleGraph &graph,
                        size_t num_paths_per_component, size_t node_window_size,
                        size_t min_node_depth, size_t max_number_of_paths_generable,
                        bool write_node_depth, std::string &node_depth,
                        const uint64_t& nthreads, const bool& ignore_paths, const bool& show_progress);

void hogwild_path_cover(handlegraph::MutablePathDeletableHandleGraph &graph,
                        double target_depth,
                        const uint64_t& nthreads, const bool& ignore_paths, const bool& show_progress);

}
}
