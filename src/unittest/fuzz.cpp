#include "catch.hpp"

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include "graph.hpp"

#include <iostream>
#include <limits>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <random>

namespace odgi {
namespace unittest {

using namespace std;
using namespace handlegraph;

TEST_CASE("Handle graphs can handle huge numbers of nodes", "[handle][packed][hashgraph]") {
    
    vector<MutablePathDeletableHandleGraph*> implementations;
    
    graph_t dg;
    implementations.push_back(&dg);

    for(MutablePathDeletableHandleGraph* implementation : implementations) {
        
        MutablePathDeletableHandleGraph& graph = *implementation;
        std::vector<handle_t> handles;
        std::mt19937 rng(87);
        for (int i = 0; i < 1000; ++i) {
            handles.push_back(graph.create_handle("A"));
        }
        std::uniform_int_distribution<std::mt19937::result_type> dist(1,1000); // distribution in range [1, 1000]
        std::uniform_int_distribution<std::mt19937::result_type> lte10(1,10);
        for (auto& handle_a : handles) {
            for (int i = 0; i < lte10(rng); ++i) {
                handle_t handle_b = graph.get_handle(dist(rng));
                graph.create_edge(handle_a, handle_b);
                REQUIRE(graph.has_edge(handle_a, handle_b));
            }
            for (int i = 0; i < lte10(rng); ++i) {
                handle_t handle_b = graph.get_handle(dist(rng));
                graph.create_edge(graph.flip(handle_a), handle_b);
                REQUIRE(graph.has_edge(graph.flip(handle_a), handle_b));
            }
        }

        for (int p = 0; p < 100; ++p) {
            //std::cerr << "path " << p << std::endl;
            path_handle_t path = graph.create_path_handle(std::to_string(p));
            handle_t last; // ...
            for (int i = 0; i < 1000; ++i) {
                //std::cerr << "occ " << i << std::endl;
                handle_t occ = handles.at(dist(rng)-1);
                graph.append_occurrence(path, occ);
                if (i > 0) {
                    graph.create_edge(last, occ);
                    REQUIRE(graph.has_edge(last, occ));
                }
                REQUIRE(graph.get_occurrence_count(path) == i+1);
                last = occ;
            }
        }
    }

}

}
}
