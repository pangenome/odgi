#pragma once

#include <unordered_map>
#include <stack>
#include <map>

#include <handlegraph/util.hpp>
#include "weakly_connected_components.hpp"

namespace odgi {
    namespace algorithms {

        using namespace handlegraph;

        constexpr size_t PATH_COVER_DEFAULT_N = 16;
        constexpr size_t PATH_COVER_DEFAULT_K = 4;

        /*
          From https://github.com/jltsiren/gbwtgraph

          Find a path cover of the graph with n paths per component and return a GBWT of the paths.
          The path cover is built greedily. Each time we extend a path, we choose the extension
          where the coverage of the k >= 2 node window is the lowest. Note that this is a maximum
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
            - If the component is not a DAG and the path is shorter than k - 1 nodes, we consider
              the coverage of individual nodes instead of windows.
            - When determining window coverage, we consider the window equivalent to its reverse
              complement.
        */
        /*gbwt::GBWT path_cover_gbwt(const HandleGraph &graph,
                                   size_t n = PATH_COVER_DEFAULT_N, size_t k = PATH_COVER_DEFAULT_K,
                                   gbwt::size_type batch_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE,
                                   gbwt::size_type sample_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL,
                                   bool show_progress = false);*/
    }

}
