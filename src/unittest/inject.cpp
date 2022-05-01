/**
 * \file
 * unittest/handle.cpp: test cases for the implementations of the HandleGraph class.
 */

#include "catch.hpp"

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include "odgi.hpp"
#include "inject.hpp"

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

TEST_CASE( "Inject works", "[inject]" ) {

    SECTION("Injection as expected") {

        graph_t g;
        string s = "CGA"; handle_t n0 = g.create_handle(s);
        s = "TTGG"; handle_t n1 = g.create_handle(s);
        s = "CCGT"; handle_t n2 = g.create_handle(s);
        s = "C"; handle_t n3 = g.create_handle(s);
        s = "GT"; handle_t n4 = g.create_handle(s);
        s = "GATAA"; handle_t n5 = g.create_handle(s);
        s = "CGG"; handle_t n6 = g.create_handle(s);
        s = "ACA"; handle_t n7 = g.create_handle(s);
        s = "GCCG"; handle_t n8 = g.create_handle(s);
        s = "ATATAAC"; handle_t n9 = g.create_handle(s);
        g.create_edge(n0,n1);
        g.create_edge(n1,n2);
        g.create_edge(n2,n3);
        g.create_edge(n2,n4);
        g.create_edge(n4,n6);
        g.create_edge(n3,n5);
        g.create_edge(n4,n5);
        g.create_edge(n5,n7);
        g.create_edge(n6,n7);
        g.create_edge(n7,n8);
        g.create_edge(n8,n9);

        path_handle_t pA = g.create_path_handle("a");
        g.append_step(pA, n0);
        g.append_step(pA, n1);
        g.append_step(pA, n2);
        g.append_step(pA, n4);
        g.append_step(pA, n6);
        g.append_step(pA, n7);
        g.append_step(pA, n8);
        g.append_step(pA, n9);

        path_handle_t pB = g.create_path_handle("b");
        g.append_step(pB, n0);
        g.append_step(pB, n1);
        g.append_step(pB, n2);
        g.append_step(pB, n3);
        g.append_step(pB, n5);
        g.append_step(pB, n7);
        g.append_step(pB, n8);
        g.append_step(pB, n9);

        path_handle_t pC = g.create_path_handle("c");
        g.append_step(pC, n0);
        g.append_step(pC, n1);
        g.append_step(pC, n2);
        g.append_step(pC, n4);
        g.append_step(pC, n5);
        g.append_step(pC, n7);
        g.append_step(pC, n8);
        g.append_step(pC, n9);

        path_handle_t pD = g.create_path_handle("d");
        g.append_step(pD, g.flip(n9));
        g.append_step(pD, g.flip(n8));
        g.append_step(pD, g.flip(n7));
        g.append_step(pD, g.flip(n5));
        g.append_step(pD, g.flip(n4));
        g.append_step(pD, g.flip(n2));
        g.append_step(pD, g.flip(n1));
        g.append_step(pD, g.flip(n0));

        ska::flat_hash_map<path_handle_t, std::vector<std::pair<interval_t, std::string>>> ivals;
        ivals[pA].push_back(std::make_pair(interval_t(0,5), "a:0-5"));
        ivals[pA].push_back(std::make_pair(interval_t(4,6), "a:4-6"));
        ivals[pA].push_back(std::make_pair(interval_t(5,10), "a:5-10"));
        ivals[pB].push_back(std::make_pair(interval_t(4,5), "b:4-5"));
        ivals[pD].push_back(std::make_pair(interval_t(4,5), "d:4-5"));
        ivals[pD].push_back(std::make_pair(interval_t(0,1), "d:0-1"));

        for (auto& i : ivals) {
            std::sort(i.second.begin(), i.second.end());
        }
        std::vector<string> ordered_ivals;
        for (auto& i : ivals) {
            for (auto& n : i.second) {
                ordered_ivals.push_back(n.second);
            }
        }

        algorithms::inject_ranges(g, ivals, ordered_ivals, true);

        //g.to_gfa(std::cout);

        REQUIRE(g.get_sequence(g.get_handle(3)) == "T");
        REQUIRE(g.get_id(g.get_handle_of_step(g.path_begin(g.get_path_handle("d:4-5")))) == 15);
        REQUIRE(g.get_sequence(g.get_handle_of_step(g.path_begin(g.get_path_handle("d:4-5")))) == "T");
        REQUIRE(g.get_step_count(g.get_path_handle("a:0-5")) == 3);
        REQUIRE(g.get_step_count(g.get_path_handle("b:4-5")) == 1);
        REQUIRE(g.get_step_count(g.get_path_handle("a:4-6")) == 2);
        REQUIRE(g.get_step_count(g.get_path_handle("a:5-10")) == 3);
        REQUIRE(g.get_step_count(g.get_path_handle("d:0-1")) == 1);
        REQUIRE(g.get_step_count(g.get_path_handle("d:4-5")) == 1);

    }
}

}
}
