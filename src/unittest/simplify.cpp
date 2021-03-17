#include "catch.hpp"

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

TEST_CASE("Graph simplification reduces a simple graph with paths", "[simplify]") {
    graph_t graph;
    handle_t n6 = graph.create_handle("TTG");
    handle_t n4 = graph.create_handle("T");
    handle_t n3 = graph.create_handle("G");
    handle_t n5 = graph.create_handle("C");
    handle_t n1 = graph.create_handle("CAAATAAG");
    handle_t n2 = graph.create_handle("A");
    graph.create_edge(n1, n2);
    graph.create_edge(n2, n3);
    graph.create_edge(n3, n4);
    graph.create_edge(n3, n5);
    graph.create_edge(n4, n5);
    graph.create_edge(n5, n6);
    path_handle_t p_x = graph.create_path_handle("x");
    path_handle_t p_y = graph.create_path_handle("y");
    for (auto& p : { p_x, p_y }) {
        for (auto& h : { n1, n2, n3, n4, n5, n6 }) {
            graph.append_step(p, h);
        }
    }
    algorithms::unchop(graph);
    // sort the graph
    graph.apply_ordering(algorithms::topological_order(&graph), true);
    graph.apply_ordering(algorithms::topological_order(&graph), true);
    graph.apply_ordering(algorithms::topological_order(&graph), true);
    SECTION("The graph is as expected") {
        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "CAAATAAGAG");
        REQUIRE(graph.get_node_count() == 3);
        step_handle_t s = graph.path_begin(p_x);
        REQUIRE(graph.get_handle_of_step(s) == graph.get_handle(1));
        s = graph.get_next_step(s);
        REQUIRE(graph.get_handle_of_step(s) == graph.get_handle(2));
        s = graph.get_next_step(s);
        REQUIRE(graph.get_handle_of_step(s) == graph.get_handle(3));
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

TEST_CASE("Unchop reduces a graph with decreasing node ids in the topological sort", "[simplify]") {
    graph_t graph;
    handle_t n1 = graph.create_handle("TTG");
    handle_t n2 = graph.create_handle("C");
    handle_t n3 = graph.create_handle("T");
    handle_t n4 = graph.create_handle("G");
    handle_t n5 = graph.create_handle("A");
    handle_t n6 = graph.create_handle("CAAATAAG");
    //handle_t n0 = graph.create_handle("TTTTTTTT");
    graph.create_edge(n2, n1);
    graph.create_edge(n3, n2);
    graph.create_edge(n4, n3);
    graph.create_edge(n5, n4);
    graph.create_edge(n6, n5);
    //graph.create_edge(n0, n4);
    algorithms::unchop(graph);
    // sort the graph
    graph.apply_ordering(algorithms::topological_order(&graph), true);
    SECTION("The graph is as expected") {
        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "CAAATAAGAGTCTTG");
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

TEST_CASE("Graph simplification reduces a graph with a self inverting -/+ loop with paths", "[simplify]") {
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
    path_handle_t p_x = graph.create_path_handle("x");
    path_handle_t p_y = graph.create_path_handle("y");
    for (auto& p : { p_x, p_y }) {
        for (auto& h : { n1, n2, n3, n4, n5, n6 }) {
            graph.append_step(p, h);
        }
    }
    path_handle_t q_y = graph.create_path_handle("q");
    for (auto& p : { q_y }) {
        for (auto& h : { n6, n5, n4, n3, n2, n1 }) {
            graph.append_step(p, graph.flip(h));
        }
    }

    algorithms::unchop(graph);

    uint64_t seen_steps = 0;
    for (auto& p : { p_x, p_y, q_y }) {
        step_handle_t begin = graph.path_begin(p);
        step_handle_t end = graph.path_end(p);
        for (step_handle_t step = begin;
             step != end;
             step = graph.get_next_step(step)) {
            handle_t h = graph.get_handle_of_step(step);
            ++seen_steps;
        }
    }

    // sort the graph
    graph.apply_ordering(algorithms::topological_order(&graph), true);

    // check that iteration still works
    uint64_t seen_steps_after_sort = 0;
    for (auto& p : { p_x, p_y, q_y }) {
        step_handle_t begin = graph.path_begin(p);
        step_handle_t end = graph.path_end(p);
        for (step_handle_t step = begin;
             step != end;
             step = graph.get_next_step(step)) {
            handle_t h = graph.get_handle_of_step(step);
            ++seen_steps_after_sort;
        }
    }

    SECTION("The graph is as expected") {
        REQUIRE(seen_steps == 6);
        REQUIRE(seen_steps == seen_steps_after_sort);
        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "CAAATAAGA");
        REQUIRE(graph.get_sequence(graph.get_handle(2)) == "GTCTTG");
        REQUIRE(graph.has_edge(graph.get_handle(1), graph.get_handle(2)));
        REQUIRE(graph.has_edge(graph.flip(graph.get_handle(2)), graph.get_handle(2)));
    }

    // sort the graph
    auto order = algorithms::topological_order(&graph);
    std::reverse(order.begin(), order.end());
    graph.apply_ordering(order, true);
    //graph.optimize(); // breaks!!

    // check that iteration still works
    uint64_t seen_steps_after_rev = 0;
    for (auto& p : { p_x, p_y, q_y }) {
        step_handle_t begin = graph.path_begin(p);
        step_handle_t end = graph.path_end(p);
        for (step_handle_t step = begin;
             step != end;
             step = graph.get_next_step(step)) {
            handle_t h = graph.get_handle_of_step(step);
            ++seen_steps_after_rev;
        }
    }

    SECTION("The graph is as expected after reverse sorting") {
        REQUIRE(seen_steps == seen_steps_after_rev);
        REQUIRE(graph.get_sequence(graph.get_handle(2)) == "CAAATAAGA");
        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "GTCTTG");
    }

}

}
}
