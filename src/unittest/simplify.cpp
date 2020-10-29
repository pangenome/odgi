#include "catch.hpp"

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include "odgi.hpp"
#include "algorithms/simple_components.hpp"
#include "algorithms/topological_sort.hpp"
#include "algorithms/unchop.hpp"

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

TEST_CASE("Graph simplification reduces a simple graph to a single node", "[simplify]") {
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
    algorithms::unchop(graph);
    // sort the graph
    graph.apply_ordering(algorithms::topological_order(&graph), true);
    SECTION("The graph is as expected") {
        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "CAAATAAGAGTCTTG");
        REQUIRE(graph.get_node_count() == 1);
    }
}

TEST_CASE("Graph simplification reduces a graph with a self loop", "[simplify]") {
    graph_t graph;
    handle_t n1 = graph.create_handle("CAAATAAG");
    handle_t n2 = graph.create_handle("A");
    handle_t n3 = graph.create_handle("G");
    handle_t n4 = graph.create_handle("T");
    handle_t n5 = graph.create_handle("C");
    handle_t n6 = graph.create_handle("TTG");
    graph.create_edge(n1, n2);
    graph.create_edge(n2, n3);
    graph.create_edge(n3, n3);
    graph.create_edge(n3, n4);
    graph.create_edge(n4, n5);
    graph.create_edge(n5, n6);
    std::vector<std::vector<handle_t>> linear_components = algorithms::simple_components(graph, 2, false, 1);
    for (auto& v : linear_components) {
        graph.combine_handles(v);
    }
    // sort the graph
    graph.apply_ordering(algorithms::topological_order(&graph), true);
    SECTION("The graph is as expected") {
        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "CAAATAAGA");
        REQUIRE(graph.get_sequence(graph.get_handle(2)) == "G");
        REQUIRE(graph.get_sequence(graph.get_handle(3)) == "TCTTG");
        REQUIRE(graph.has_edge(graph.get_handle(1), graph.get_handle(2)));
        REQUIRE(graph.has_edge(graph.get_handle(2), graph.get_handle(2)));
        REQUIRE(graph.has_edge(graph.get_handle(2), graph.get_handle(3)));
    }
}

TEST_CASE("Unchop reduces a graph with a self inverting +/- loop", "[simplify]") {
    graph_t graph;
    handle_t n1 = graph.create_handle("CAAATAAG");
    handle_t n2 = graph.create_handle("A");
    handle_t n3 = graph.create_handle("G");
    handle_t n4 = graph.create_handle("T");
    handle_t n5 = graph.create_handle("C");
    handle_t n6 = graph.create_handle("TTG");
    graph.create_edge(n1, n2);
    graph.create_edge(n2, n3);
    graph.create_edge(n3, graph.flip(n3));
    graph.create_edge(n3, n4);
    graph.create_edge(n4, n5);
    graph.create_edge(n5, n6);
    algorithms::unchop(graph);
    // sort the graph
    graph.apply_ordering(algorithms::topological_order(&graph), true);
    SECTION("The graph is as expected") {
        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "CAAATAAGAG");
        REQUIRE(graph.get_sequence(graph.get_handle(2)) == "TCTTG");
        REQUIRE(graph.has_edge(graph.get_handle(1), graph.get_handle(2)));
        REQUIRE(graph.has_edge(graph.get_handle(1), graph.flip(graph.get_handle(1))));
    }
}

TEST_CASE("Unchop reduces a graph that makes a simple loop", "[simplify]") {
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
    graph.create_edge(n6, n1);
    algorithms::unchop(graph);
    // sort the graph
    graph.apply_ordering(algorithms::topological_order(&graph), true);
    SECTION("The graph is as expected") {
        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "CAAATAAGAGTCTTG");
        REQUIRE(graph.has_edge(graph.get_handle(1), graph.get_handle(1)));
    }
}

TEST_CASE("Graph simplification reduces a graph with a self inverting -/+ loop", "[simplify]") {
    graph_t graph;
    handle_t n1 = graph.create_handle("CAAATAAG");
    handle_t n2 = graph.create_handle("A");
    handle_t n3 = graph.create_handle("G");
    handle_t n4 = graph.create_handle("T");
    handle_t n5 = graph.create_handle("C");
    handle_t n6 = graph.create_handle("TTG");
    graph.create_edge(n1, n2);
    graph.create_edge(n2, n3);
    graph.create_edge(graph.flip(n3), n3);
    graph.create_edge(n3, n4);
    graph.create_edge(n4, n5);
    graph.create_edge(n5, n6);
    algorithms::unchop(graph);
    // sort the graph
    graph.apply_ordering(algorithms::topological_order(&graph), true);
    SECTION("The graph is as expected") {
        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "CAAATAAGA");
        REQUIRE(graph.get_sequence(graph.get_handle(2)) == "GTCTTG");
        REQUIRE(graph.has_edge(graph.get_handle(1), graph.get_handle(2)));
        REQUIRE(graph.has_edge(graph.flip(graph.get_handle(2)), graph.get_handle(2)));
    }
}

}
}
