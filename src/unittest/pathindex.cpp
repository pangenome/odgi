/**
 * \file
 * unittest/pathindex.cpp: test cases for the implementations of the XP class.
 */

#include "catch.hpp"

#include <handlegraph/util.hpp>
#include "odgi.hpp"
#include "algorithms/xp.hpp"

namespace odgi {
    namespace unittest {

        using namespace std;
        using namespace handlegraph;
        using namespace xp;

        TEST_CASE("XP construction.", "[path_index]") {
            // the name of the path in the graph
            string f = "5";
            graph_t graph;
            handle_t n1 = graph.create_handle("AGGA");
            handle_t n2 = graph.create_handle("A");
            handle_t n3 = graph.create_handle("TC");
            handle_t n4 = graph.create_handle("TCTCAGG");
            graph.create_edge(n1, n2);
            graph.create_edge(n2, n3);
            graph.create_edge(n2, n4);
            graph.create_edge(n3, n4);

            graph.create_path_handle(f, false);
            path_handle_t five = graph.get_path_handle(f);
            graph.append_step(five, n1);
            graph.append_step(five, n3);
            graph.append_step(five, n4);

            // just test if the graph generation worked roughly
            SECTION("The graph is as expected") {
                REQUIRE(graph.get_node_count() == 4);
                REQUIRE(graph.has_path(f));
                REQUIRE(graph.get_length(n1) == 4);
            }

            XP path_index;
            path_index.from_handle_graph(graph);

            SECTION("The index mirrors the actual graph") {
                REQUIRE(path_index.path_count == graph.get_path_count());
                REQUIRE(path_index.has_path(f));

                // TODO
                /*
                path_handle_t p_h_i = path_index.get_path_handle(f);
                step_handle_t s_h_i = path_index.get_step_at_position(p_h_i, 12);
                step_handle_t s_h_g = graph. // TODO @ekg I don't know how to do this.
                REQUIRE(path_index.get_path_handle_of_step(s_h_i) == graph.get_path_handle_of_step(s_h_g));
                */

                // FIXME @ekg This does not run through.
                // REQUIRE(as_integer(path_index.get_path_handle("5")) == as_integer(graph.get_path_handle("5")));
            }

            SECTION("The index has path and position") {
                REQUIRE(!path_index.has_path("4"));
                REQUIRE(path_index.has_position("5", 5));
                REQUIRE(path_index.has_position("5", 1));
                REQUIRE(path_index.has_position("5", 13));
                REQUIRE(!path_index.has_position("5", 45));
                REQUIRE(!path_index.has_position("4", 5));
                REQUIRE(!path_index.has_position("4", 45));
            }

        }
    }
}