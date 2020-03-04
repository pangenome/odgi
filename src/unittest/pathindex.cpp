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

        TEST_CASE("XP construction, serialization and loading.", "[pathindex]") {

            graph_t graph;
            handle_t n1 = graph.create_handle("AGGA");
            handle_t n2 = graph.create_handle("A");
            handle_t n3 = graph.create_handle("TC");
            handle_t n3_m = graph.flip(n3);
            handle_t n4 = graph.create_handle("TCTCAGG");
            graph.create_edge(n1, n2);
            graph.create_edge(n2, n3);
            graph.create_edge(n2, n4);
            graph.create_edge(n3, n4);
            graph.create_edge(n1, n4);
            graph.create_edge(n4, n3_m);
            // TODO REMOVE THIS
            graph.create_edge(n4, n3);

            graph.create_path_handle("5", false);
            path_handle_t five = graph.get_path_handle("5");
            graph.append_step(five, n1);
            graph.append_step(five, n3);
            graph.append_step(five, n4);

            graph.create_path_handle("5-", false);
            path_handle_t five_m = graph.get_path_handle("5-");
            graph.append_step(five_m, n1);
            graph.append_step(five_m, n4);
            graph.append_step(five_m, n3);

            // just test if the graph generation worked roughly
            SECTION("The graph is as expected") {
                REQUIRE(graph.get_node_count() == 4);
                REQUIRE(graph.has_path("5"));
                REQUIRE(graph.has_path("5-"));
                REQUIRE(!graph.has_path("5+"));
                REQUIRE(graph.get_length(n1) == 4);
                REQUIRE(graph.get_length(n3_m) == 2);
                REQUIRE(graph.get_length(n3) == 2);
                REQUIRE(graph.get_sequence(n1) == "AGGA");
                REQUIRE(graph.get_sequence(n3_m) == "GA");
            }

            XP path_index;
            path_index.from_handle_graph(graph);

            SECTION("The index mirrors the actual graph") {
                REQUIRE(path_index.path_count == graph.get_path_count());
                REQUIRE(path_index.has_path("5"));
                REQUIRE(path_index.has_path("5-"));
                REQUIRE(!path_index.has_path("5+"));

                // TODO
                /*
                path_handle_t p_h_i = path_index.get_path_handle(f);
                step_handle_t s_h_i = path_index.get_step_at_position(p_h_i, 12);
                step_handle_t s_h_g = graph. // TODO @ekg I don't know how to do this.
                REQUIRE(path_index.get_path_handle_of_step(s_h_i) == graph.get_path_handle_of_step(s_h_g));
                */

                // We have to add +1 in the graph space, because of the different handle implementation in XP
                // Explanation by Erik: https://github.com/vgteam/odgi/pull/82#discussion_r387129361.
                REQUIRE(as_integer(path_index.get_path_handle("5")) == as_integer(graph.get_path_handle("5"))+ 1);
                REQUIRE(as_integer(path_index.get_path_handle("5-")) == as_integer(graph.get_path_handle("5-"))+ 1);

                SECTION("The XPPaths mirror the actual paths in the graph") {

                }
            }

            SECTION("The index has path and position") {
                REQUIRE(!path_index.has_path("4"));
                REQUIRE(path_index.has_position("5", 5));
                REQUIRE(path_index.has_position("5", 1));
                REQUIRE(path_index.has_position("5", 13));
                REQUIRE(!path_index.has_position("5", 45));
                REQUIRE(!path_index.has_position("4", 5));
                REQUIRE(!path_index.has_position("4", 45));
                REQUIRE(path_index.has_position("5-", 5));
                REQUIRE(path_index.has_position("5-", 1));
                REQUIRE(path_index.has_position("5-", 13));
                REQUIRE(!path_index.has_position("5-", 45));
                REQUIRE(!path_index.has_position("4-", 5));
                REQUIRE(!path_index.has_position("4+", 45));
            }

            SECTION("Retrieving pangenome position from constructed index") {
                REQUIRE(path_index.get_pangenome_pos("5", 1) == 1);
                REQUIRE(path_index.get_pangenome_pos("5", 2) == 2);
                REQUIRE(path_index.get_pangenome_pos("5", 13) == 14);
                REQUIRE(path_index.get_pangenome_pos("5", 5) == 6);
                REQUIRE(path_index.get_pangenome_pos("5", 12) == 13);
                REQUIRE(path_index.get_pangenome_pos("5", 24) == 0);
                REQUIRE(path_index.get_pangenome_pos("4", 1) == 0);
                REQUIRE(path_index.get_pangenome_pos("5-", 1) == 1);
                REQUIRE(path_index.get_pangenome_pos("5-", 2) == 2);
                REQUIRE(path_index.get_pangenome_pos("5-", 6) == 9);
                REQUIRE(path_index.get_pangenome_pos("5-", 13) == 14);
                REQUIRE(path_index.get_pangenome_pos("5-", 5) == 6);
                REQUIRE(path_index.get_pangenome_pos("5-", 12) == 13);
                REQUIRE(path_index.get_pangenome_pos("5-", 24) == 0);
                REQUIRE(path_index.get_pangenome_pos("4", 1) == 0);
            }

            // Write index to temporary file in preparation for the next test section.
            std::string basename = temp_file::create();
            std::ofstream out;
            out.open(basename + "unittest_pathindex.xp");
            path_index.serialize_members(out);
            out.close();

            XP loaded_path_index;
            std::ifstream in;
            in.open(basename + "unittest_pathindex.xp");
            loaded_path_index.load(in);
            in.close();

            SECTION("The loaded index mirrors the actual graph") {
                REQUIRE(loaded_path_index.path_count == graph.get_path_count());
                REQUIRE(loaded_path_index.has_path("5"));

                // TODO
                /*
                path_handle_t p_h_i = loaded_path_index.get_path_handle(f);
                step_handle_t s_h_i = loaded_path_index.get_step_at_position(p_h_i, 12);
                step_handle_t s_h_g = graph. // TODO @ekg I don't know how to do this.
                REQUIRE(loaded_path_index.get_path_handle_of_step(s_h_i) == graph.get_path_handle_of_step(s_h_g));
                */

                // We have to add +1 in the graph space, because of the different handle implementation in XP
                // Explanation by Erik: https://github.com/vgteam/odgi/pull/82#discussion_r387129361.
                REQUIRE(as_integer(loaded_path_index.get_path_handle("5")) == as_integer(graph.get_path_handle("5"))+ 1);
            }

            SECTION("The loaded index mirrors the actual index") {
                REQUIRE(loaded_path_index.path_count == graph.get_path_count());
                REQUIRE(loaded_path_index.has_path("5"));
                REQUIRE(as_integer(loaded_path_index.get_path_handle("5")) == as_integer(path_index.get_path_handle("5")));
                REQUIRE(loaded_path_index.get_path_handle("5") == path_index.get_path_handle("5"));

                REQUIRE(loaded_path_index.get_path("5").min_handle == path_index.get_path("5").min_handle);
            }

            SECTION("The loaded index has path and position") {
                REQUIRE(!loaded_path_index.has_path("4"));
                REQUIRE(loaded_path_index.has_position("5", 5));
                REQUIRE(loaded_path_index.has_position("5", 1));
                REQUIRE(loaded_path_index.has_position("5", 13));
                REQUIRE(!loaded_path_index.has_position("5", 45));
                REQUIRE(!loaded_path_index.has_position("4", 5));
                REQUIRE(!loaded_path_index.has_position("4", 45));
            }

            SECTION("Retrieving pangenome position from loaded index") {
                REQUIRE(loaded_path_index.get_pangenome_pos("5", 1) == 1);
                REQUIRE(loaded_path_index.get_pangenome_pos("5", 2) == 2);
                REQUIRE(loaded_path_index.get_pangenome_pos("5", 13) == 14);
                REQUIRE(loaded_path_index.get_pangenome_pos("5", 5) == 6);
                REQUIRE(loaded_path_index.get_pangenome_pos("5", 12) == 13);
                REQUIRE(loaded_path_index.get_pangenome_pos("5", 24) == 0);
                REQUIRE(loaded_path_index.get_pangenome_pos("4", 1) == 0);
            }
        }
    }
}