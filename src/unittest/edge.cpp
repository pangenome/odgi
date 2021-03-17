/**
 * \file
 * unittest/handle.cpp: test cases for the implementations of the HandleGraph class.
 */

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

TEST_CASE( "Follow edges utility works", "[edges]" ) {

    SECTION("From node to left node") {

        SECTION("Forward to forward") {
            graph_t graph;
            handle_t node = graph.create_handle("A");
            handle_t left_node = graph.create_handle("T");
            handle_t right_node = graph.create_handle("C");
            graph.create_edge(left_node, node);
            vector<handle_t> handles;
            // follow all edges that can go to the left node
            graph.follow_edges(node, true, [&](const handle_t& h) {
                //std::cerr << "node_left rank: " << number_bool_packing::unpack_number(h) << ", node_left id: "
                //<< graph.get_id(h) << std::endl;
                handles.push_back(h);
            });
            REQUIRE(1 == handles.size());
            //std::cerr << "node rank: " << number_bool_packing::unpack_number(node) << ", node id: "
            //    << graph.get_id(node) << std::endl;
        }

        SECTION("Forward to backward") {
            graph_t graph;
            handle_t node = graph.create_handle("A");
            handle_t left_node = graph.create_handle("T");
            handle_t right_node = graph.create_handle("C");
            graph.create_edge(left_node, graph.flip(node));
            vector<handle_t> handles;
            graph.follow_edges(node, false, [&](const handle_t& h) {
                //std::cerr << "node_left rank: " << number_bool_packing::unpack_number(h) << ", node_left id: "
                //<< graph.get_id(h) << std::endl;
                handles.push_back(h);
            });
            REQUIRE(1 == handles.size());
        }
        SECTION("Backward to forward") {
            graph_t graph;
            handle_t node = graph.create_handle("A");
            handle_t left_node = graph.create_handle("T");
            handle_t right_node = graph.create_handle("C");
            graph.create_edge(graph.flip(left_node), node);
            vector<handle_t> handles;
            graph.follow_edges(node, true, [&](const handle_t& h) {
                //std::cerr << "node_left rank: " << number_bool_packing::unpack_number(h) << ", node_left id: "
                //<< graph.get_id(h) << std::endl;
                handles.push_back(h);
            });
            REQUIRE(1 == handles.size());
        }
        SECTION("Backward to backward") {
            graph_t graph;
            handle_t node = graph.create_handle("A");
            handle_t left_node = graph.create_handle("T");
            handle_t right_node = graph.create_handle("C");
            graph.create_edge(graph.flip(left_node), graph.flip(node));
            vector<handle_t> handles;
            graph.follow_edges(node, false, [&](const handle_t& h) {
                //std::cerr << "node_left rank: " << number_bool_packing::unpack_number(h) << ", node_left id: "
                //<< graph.get_id(h) << std::endl;
                handles.push_back(h);
            });
            REQUIRE(1 == handles.size());
        }
    }

    SECTION("From node to right node") {

        SECTION("Forward to forward") {
            graph_t graph;
            handle_t node = graph.create_handle("A");
            handle_t left_node = graph.create_handle("T");
            handle_t right_node = graph.create_handle("C");
            graph.create_edge(node, right_node);
            vector<handle_t> handles;
            // follow all edges that can go to the right node
            graph.follow_edges(node, false, [&](const handle_t& h) {
                //std::cerr << "node_left rank: " << number_bool_packing::unpack_number(h) << ", node_left id: "
                //<< graph.get_id(h) << std::endl;
                handles.push_back(h);
            });
            REQUIRE(1 == handles.size());
            //std::cerr << "node rank: " << number_bool_packing::unpack_number(node) << ", node id: "
            //    << graph.get_id(node) << std::endl;
        }

        SECTION("Forward to backward") {
            graph_t graph;
            handle_t node = graph.create_handle("A");
            handle_t left_node = graph.create_handle("T");
            handle_t right_node = graph.create_handle("C");
            graph.create_edge(node, graph.flip(right_node));
            vector<handle_t> handles;
            graph.follow_edges(node, false, [&](const handle_t& h) {
                //std::cerr << "node_left rank: " << number_bool_packing::unpack_number(h) << ", node_left id: "
                //<< graph.get_id(h) << std::endl;
                handles.push_back(h);
            });
            REQUIRE(1 == handles.size());
        }
        SECTION("Backward to forward") {
            graph_t graph;
            handle_t node = graph.create_handle("A");
            handle_t left_node = graph.create_handle("T");
            handle_t right_node = graph.create_handle("C");
            graph.create_edge(graph.flip(node), right_node);
            vector<handle_t> handles;
            graph.follow_edges(node, true, [&](const handle_t& h) {
                //std::cerr << "node_left rank: " << number_bool_packing::unpack_number(h) << ", node_left id: "
                //<< graph.get_id(h) << std::endl;
                handles.push_back(h);
            });
            REQUIRE(1 == handles.size());
        }
        SECTION("Backward to backward") {
            graph_t graph;
            handle_t node = graph.create_handle("A");
            handle_t left_node = graph.create_handle("T");
            handle_t right_node = graph.create_handle("C");
            graph.create_edge(graph.flip(node), graph.flip(right_node));
            vector<handle_t> handles;
            graph.follow_edges(node, true, [&](const handle_t& h) {
                //std::cerr << "node_left rank: " << number_bool_packing::unpack_number(h) << ", node_left id: "
                //<< graph.get_id(h) << std::endl;
                handles.push_back(h);
            });
            REQUIRE(1 == handles.size());
        }
    }

    SECTION("From node to node") {

        SECTION("Forward to forward") {
            graph_t graph;
            handle_t node = graph.create_handle("A");
            handle_t left_node = graph.create_handle("T");
            handle_t right_node = graph.create_handle("C");
            graph.create_edge(node, node);
            vector<handle_t> handles;
            // follow all edges that can go to from node to node
            graph.follow_edges(node, false, [&](const handle_t& h) {
                handles.push_back(h);
            });
            REQUIRE(1 == handles.size());

            graph.follow_edges(node, true, [&](const handle_t& h) {

                handles.push_back(h);
            });
            REQUIRE(2 == handles.size());
        }

        SECTION("Forward to backward") {
            graph_t graph;
            handle_t node = graph.create_handle("A");
            handle_t left_node = graph.create_handle("T");
            handle_t right_node = graph.create_handle("C");
            graph.create_edge(node, graph.flip(node));
            vector<handle_t> handles;
            graph.follow_edges(node, false, [&](const handle_t& h) {
                //std::cerr << "node_left rank: " << number_bool_packing::unpack_number(h) << ", node_left id: "
                //<< graph.get_id(h) << std::endl;
                handles.push_back(h);
            });
            REQUIRE(1 == handles.size());
            graph.follow_edges(node, true, [&](const handle_t& h) {
                //std::cerr << "node_left rank: " << number_bool_packing::unpack_number(h) << ", node_left id: "
                //<< graph.get_id(h) << std::endl;
                handles.push_back(h);
            });
            REQUIRE(1 == handles.size());
        }
        SECTION("Backward to forward") {
            graph_t graph;
            handle_t node = graph.create_handle("A");
            handle_t left_node = graph.create_handle("T");
            handle_t right_node = graph.create_handle("C");
            graph.create_edge(graph.flip(node), node);
            vector<handle_t> handles;
            graph.follow_edges(node, true, [&](const handle_t& h) {
                //std::cerr << "node_left rank: " << number_bool_packing::unpack_number(h) << ", node_left id: "
                //<< graph.get_id(h) << std::endl;
                handles.push_back(h);
            });
            REQUIRE(1 == handles.size());
            graph.follow_edges(node, false, [&](const handle_t& h) {
                handles.push_back(h);
            });
            REQUIRE(1 == handles.size());
        }
        SECTION("Backward to backward") {
            graph_t graph;
            handle_t node = graph.create_handle("A");
            handle_t left_node = graph.create_handle("T");
            handle_t right_node = graph.create_handle("C");
            graph.create_edge(graph.flip(node), graph.flip(node));
            vector<handle_t> handles;
            graph.follow_edges(node, false, [&](const handle_t& h) {
                handles.push_back(h);
            });
            REQUIRE(1 == handles.size());
            graph.follow_edges(node, true, [&](const handle_t& h) {
                handles.push_back(h);
            });
            REQUIRE(2 == handles.size());
        }
    }
}

}
}