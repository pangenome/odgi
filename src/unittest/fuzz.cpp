#include "catch.hpp"

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include "odgi.hpp"

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

TEST_CASE("Large random handle with high path depth work", "[fuzz]") {

    vector<MutablePathDeletableHandleGraph*> implementations;

    graph_t dg;
    implementations.push_back(&dg);

    for(MutablePathDeletableHandleGraph* implementation : implementations) {

        MutablePathDeletableHandleGraph& graph = *implementation;
        std::vector<handle_t> handles;
        std::mt19937 rng(87);
        int num_handles = 1e6;
        std::uniform_int_distribution<std::mt19937::result_type> dist(1,num_handles);
        std::uniform_int_distribution<std::mt19937::result_type> lte5(1,5);
        std::uniform_int_distribution<std::mt19937::result_type> lte1000(1,1000);
        for (int i = 0; i < num_handles; ++i) {
            std::string s('A', lte1000(rng));
            handles.push_back(graph.create_handle(s));
            REQUIRE(graph.get_sequence(handles.back()) == s);
        }
        for (auto& handle_a : handles) {
            for (int i = 0; i < lte5(rng); ++i) {
                handle_t handle_b = graph.get_handle(dist(rng));
                graph.create_edge(handle_a, handle_b);
                REQUIRE(graph.has_edge(handle_a, handle_b));
            }
            for (int i = 0; i < lte5(rng); ++i) {
                handle_t handle_b = graph.get_handle(dist(rng));
                graph.create_edge(graph.flip(handle_a), handle_b);
                REQUIRE(graph.has_edge(graph.flip(handle_a), handle_b));
            }
        }

        for (int p = 0; p < 100; ++p) {
            path_handle_t path = graph.create_path_handle(std::to_string(p));
            handle_t last; // ...
            for (int i = 0; i < 100000; ++i) {
                handle_t occ = handles.at(dist(rng)-1);
                graph.append_step(path, occ);
                if (i > 0) {
                    graph.create_edge(last, occ);
                    REQUIRE(graph.has_edge(last, occ));
                }
                REQUIRE(graph.get_step_count(path) == i+1);
                last = occ;
            }
        }
    }

}

}
}
