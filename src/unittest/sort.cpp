#include "catch.hpp"

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include "odgi.hpp"
#include "algorithms/topological_sort.hpp"

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

TEST_CASE("Sorting a graph without paths", "[sort]") {
    graph_t graph;
    handle_t n1 = graph.create_handle("CAAATAAG");
    handle_t n2 = graph.create_handle("A");
    handle_t n3 = graph.create_handle("G");
    handle_t n4 = graph.create_handle("T");
    handle_t n5 = graph.create_handle("C");
    handle_t n6 = graph.create_handle("TTG");
    graph.create_edge(n1, n2);
    graph.create_edge(n2, n3);
    graph.create_edge(n3, n4);
    graph.create_edge(n4, n5);
    graph.create_edge(n5, n6);
    // sort the graph
    auto order = algorithms::topological_order(&graph);
    graph.apply_ordering(order, true);
    SECTION("The graph is as expected when sorted") {
        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "CAAATAAG");
        REQUIRE(graph.get_sequence(graph.get_handle(6)) == "TTG");
        REQUIRE(graph.get_node_count() == 6);
    }
    std::reverse(order.begin(), order.end());
    graph.apply_ordering(order, true);
    SECTION("The graph is as expected when reversed") {
        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "TTG");
        REQUIRE(graph.get_sequence(graph.get_handle(6)) == "CAAATAAG");
        REQUIRE(graph.get_node_count() == 6);
    }
}

TEST_CASE("Sorting a graph with paths", "[sort]") {
    graph_t graph;
    handle_t n1 = graph.create_handle("CAAATAAG");
    handle_t n2 = graph.create_handle("A");
    handle_t n3 = graph.create_handle("G");
    handle_t n4 = graph.create_handle("T");
    handle_t n5 = graph.create_handle("C");
    handle_t n6 = graph.create_handle("TTG");
    graph.create_edge(n1, n2);
    graph.create_edge(n2, n3);
    graph.create_edge(n3, n4);
    graph.create_edge(n4, n5);
    graph.create_edge(n5, n6);
    std::vector<path_handle_t> paths;
    for (uint64_t i = 0; i < 10; ++i) {
        auto path = graph.create_path_handle("x" + std::to_string(i));
        graph.append_step(path, n1);
        graph.append_step(path, n2);
        graph.append_step(path, n3);
        graph.append_step(path, n4);
        graph.append_step(path, n5);
        graph.append_step(path, n6);
        paths.push_back(path);
    }

    auto test_path =
        [&](const path_handle_t& p) {
            auto& path_meta = graph.path_metadata(p);
            uint64_t i = 0;
            graph.for_each_step_in_path(p, [&](const step_handle_t& step) {
                                               handle_t h = graph.get_handle_of_step(step);
                                               ++i;
                                           });
            REQUIRE(i == path_meta.length);
        };

    graph.for_each_path_handle(test_path);

    // sort the graph
    auto order = algorithms::topological_order(&graph);
    graph.apply_ordering(order, true);
    graph.for_each_path_handle(test_path);

    graph.for_each_path_handle(
        [&](const path_handle_t& p) {
            auto& path_meta = graph.path_metadata(p);
            uint64_t i = 0;
            graph.for_each_step_in_path(p, [&](const step_handle_t& step) {
                                         handle_t h = graph.get_handle_of_step(step);
                                         ++i;
                                     });
            REQUIRE(i == path_meta.length);
        });
    
    SECTION("The graph is as expected when sorted") {
        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "CAAATAAG");
        REQUIRE(graph.get_sequence(graph.get_handle(6)) == "TTG");
        REQUIRE(graph.get_node_count() == 6);
        for (auto& path : paths) {
            REQUIRE(graph.get_handle_of_step(graph.path_begin(path)) == n1);
            std::string x;
            graph.for_each_step_in_path(
                path,
                [&](const step_handle_t& step) {
                    x.append(graph.get_sequence(graph.get_handle_of_step(step)));
                });
            REQUIRE(x == "CAAATAAGAGTCTTG");
        }
    }
    std::reverse(order.begin(), order.end());
    graph.apply_ordering(order, true);

    graph.for_each_path_handle(test_path);

    SECTION("The graph is as expected when reversed") {
        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "TTG");
        REQUIRE(graph.get_sequence(graph.get_handle(6)) == "CAAATAAG");
        REQUIRE(graph.get_node_count() == 6);
    }
}

}
}
