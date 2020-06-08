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

            graph.create_path_handle("5-m", false);
            path_handle_t five_m_m = graph.get_path_handle("5-m");
            graph.append_step(five_m_m, n1);
            graph.append_step(five_m_m, n4);
            graph.append_step(five_m_m, n3_m);

            // just test if the graph generation worked roughly
            SECTION("The graph is as expected") {
                REQUIRE(graph.get_node_count() == 4);
                REQUIRE(graph.has_path("5"));
                REQUIRE(graph.has_path("5-"));
                REQUIRE(!graph.has_path("5+"));
                REQUIRE(graph.has_path("5-m"));
                REQUIRE(graph.get_length(n1) == 4);
                REQUIRE(graph.get_length(n3_m) == 2);
                REQUIRE(graph.get_length(n3) == 2);
                REQUIRE(graph.get_sequence(n1) == "AGGA");
                REQUIRE(graph.get_sequence(n3_m) == "GA");
            }

            XP path_index;
            // graph.destroy_handle(n2); // this leads to a 'exit(1)' as we don't create an index of an optimized graph anymore
            path_index.from_handle_graph(graph);

            SECTION("The index mirrors the actual graph") {
                REQUIRE(path_index.path_count == graph.get_path_count());
                REQUIRE(path_index.has_path("5"));
                REQUIRE(path_index.has_path("5-"));
                REQUIRE(path_index.has_path("5-m"));
                REQUIRE(!path_index.has_path("5+"));

                REQUIRE(as_integer(path_index.get_path_handle("5")) == as_integer(graph.get_path_handle("5")));
                REQUIRE(as_integer(path_index.get_path_handle("5-")) == as_integer(graph.get_path_handle("5-")));
                REQUIRE(as_integer(path_index.get_path_handle("5-m")) == as_integer(graph.get_path_handle("5-m")));
            }

            SECTION("The index has path and position") {
                REQUIRE(!path_index.has_path("4"));
                REQUIRE(path_index.has_position("5", 4));
                REQUIRE(path_index.has_position("5", 0));
                REQUIRE(path_index.has_position("5", 12));
                REQUIRE(!path_index.has_position("5", 44));
                REQUIRE(!path_index.has_position("4", 4));
                REQUIRE(!path_index.has_position("4", 44));
                REQUIRE(path_index.has_position("5-", 4));
                REQUIRE(path_index.has_position("5-", 0));
                REQUIRE(path_index.has_position("5-", 12));
                REQUIRE(!path_index.has_position("5-", 44));
                REQUIRE(!path_index.has_position("4-", 4));
                REQUIRE(!path_index.has_position("4+", 44));
                REQUIRE(path_index.has_position("5-m", 4));
                REQUIRE(path_index.has_position("5-m", 0));
                REQUIRE(path_index.has_position("5-m", 12));
                REQUIRE(!path_index.has_position("5-m", 44));
            }

            SECTION("Retrieving pangenome position from constructed index") {
                REQUIRE(path_index.get_pangenome_pos("5", 0) == 0);
                REQUIRE(path_index.get_pangenome_pos("5", 1) == 1);
                REQUIRE(path_index.get_pangenome_pos("5", 12) == 13);
                REQUIRE(path_index.get_pangenome_pos("5", 4) == 5);
                REQUIRE(path_index.get_pangenome_pos("5", 11) == 12);
                REQUIRE(path_index.get_pangenome_pos("5-", 0) == 0);
                REQUIRE(path_index.get_pangenome_pos("5-", 1) == 1);
                REQUIRE(path_index.get_pangenome_pos("5-", 5) == 8);
                REQUIRE(path_index.get_pangenome_pos("5-", 12) == 6);
                REQUIRE(path_index.get_pangenome_pos("5-", 4) == 7);
                REQUIRE(path_index.get_pangenome_pos("5-", 11) == 5);
                REQUIRE(path_index.get_pangenome_pos("5-m", 0) == 0);
                REQUIRE(path_index.get_pangenome_pos("5-m", 1) == 1);
                REQUIRE(path_index.get_pangenome_pos("5-m", 5) == 8);
                REQUIRE(path_index.get_pangenome_pos("5-m", 12) == 5);
                REQUIRE(path_index.get_pangenome_pos("5-m", 4) == 7);
                REQUIRE(path_index.get_pangenome_pos("5-m", 11) == 6);

                /// The program will exit(1) here.
                // REQUIRE(path_index.get_pangenome_pos("5-", 24) == 0);
                // REQUIRE(path_index.get_pangenome_pos("4", 1) == 0);
                // REQUIRE(path_index.get_pangenome_pos("5", 24) == 0);
                // REQUIRE(path_index.get_pangenome_pos("4", 1) == 0);
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
                REQUIRE(loaded_path_index.has_path("5-"));
                REQUIRE(loaded_path_index.has_path("5-m"));
                REQUIRE(!loaded_path_index.has_path("5+"));

                REQUIRE(as_integer(loaded_path_index.get_path_handle("5")) == as_integer(graph.get_path_handle("5")));
                REQUIRE(as_integer(loaded_path_index.get_path_handle("5-")) == as_integer(graph.get_path_handle("5-")));
                REQUIRE(as_integer(loaded_path_index.get_path_handle("5-m")) == as_integer(graph.get_path_handle("5-m")));
            }

            SECTION("The loaded index mirrors the actual index") {
                REQUIRE(loaded_path_index.path_count == graph.get_path_count());
                REQUIRE(loaded_path_index.has_path("5"));
                REQUIRE(as_integer(loaded_path_index.get_path_handle("5")) == as_integer(path_index.get_path_handle("5")));
                REQUIRE(loaded_path_index.get_path_handle("5") == path_index.get_path_handle("5"));
                REQUIRE(loaded_path_index.get_path("5").min_handle == path_index.get_path("5").min_handle);
                // compare pn_iv XP
                sdsl::int_vector<> pn_iv_l = loaded_path_index.get_pn_iv();
                sdsl::int_vector<> pn_iv = path_index.get_pn_iv();
                for (size_t i = 0; i < loaded_path_index.get_pn_iv().size(); i++) {
                    REQUIRE(pn_iv_l[i] == pn_iv[i]);
                }
                // compare pos_map_iv XP
                sdsl::enc_vector<> pos_map_iv_l = loaded_path_index.get_pos_map_iv();
                sdsl::enc_vector<> pos_map_iv = path_index.get_pos_map_iv();
                for (size_t i = 0; i < loaded_path_index.get_pos_map_iv().size(); i++) {
                    REQUIRE(pos_map_iv_l[i] == pos_map_iv[i]);
                }

                // compare all fields of all XPPATHs
                std::vector<XPPath*> xppaths_lpi = loaded_path_index.get_paths();
                std::vector<XPPath*> xppaths_pi = loaded_path_index.get_paths();
                for (size_t i = 0; i < xppaths_pi.size(); i++) {
                    XPPath* xppath_lpi = xppaths_lpi[i];
                    XPPath* xppath_pi = xppaths_pi[i];
                    REQUIRE(xppath_lpi->is_circular == xppath_pi->is_circular);
                    REQUIRE(xppath_lpi->min_handle == xppath_pi->min_handle);
                    REQUIRE(xppath_lpi->handles.size() == xppath_pi->handles.size());
                    REQUIRE(xppath_lpi->offsets.size() == xppath_pi->offsets.size());
                    for (size_t j = 0; j < xppath_lpi->handles.size(); j++) {
                        REQUIRE(xppath_lpi->handles[j] == xppath_pi->handles[j]);
                    }
                    for (size_t j = 0; j < xppath_lpi->offsets.size(); j++) {
                        REQUIRE(xppath_lpi->offsets[j] == xppath_pi->offsets[j]);
                    }
                    for (size_t j = 0; j < xppath_lpi->offsets_rank.size(); j++) {
                        REQUIRE(xppath_lpi->offsets_rank(j) == xppath_pi->offsets_rank(j));
                    }
                    for (size_t j = 0; j < xppath_lpi->offsets_select.size(); j++) {
                        REQUIRE(xppath_lpi->offsets_select(j) == xppath_pi->offsets_select(j));
                    }
                }
            }

            SECTION("The loaded index has path and position") {
                REQUIRE(!loaded_path_index.has_path("4"));
                REQUIRE(loaded_path_index.has_position("5", 4));
                REQUIRE(loaded_path_index.has_position("5", 0));
                REQUIRE(loaded_path_index.has_position("5", 12));
                REQUIRE(!loaded_path_index.has_position("5", 44));
                REQUIRE(!loaded_path_index.has_position("4", 4));
                REQUIRE(!loaded_path_index.has_position("4", 44));
                REQUIRE(loaded_path_index.has_position("5-", 4));
                REQUIRE(loaded_path_index.has_position("5-", 0));
                REQUIRE(loaded_path_index.has_position("5-", 12));
                REQUIRE(!loaded_path_index.has_position("5-", 44));
                REQUIRE(!loaded_path_index.has_position("4-", 4));
                REQUIRE(!loaded_path_index.has_position("4+", 44));
                REQUIRE(loaded_path_index.has_position("5-m", 4));
                REQUIRE(loaded_path_index.has_position("5-m", 0));
                REQUIRE(loaded_path_index.has_position("5-m", 12));
                REQUIRE(!loaded_path_index.has_position("5-m", 44));
            }

            SECTION("Retrieving pangenome position from loaded index") {
                REQUIRE(loaded_path_index.get_pangenome_pos("5", 0) == 0);
                REQUIRE(loaded_path_index.get_pangenome_pos("5", 1) == 1);
                REQUIRE(loaded_path_index.get_pangenome_pos("5", 12) == 13);
                REQUIRE(loaded_path_index.get_pangenome_pos("5", 4) == 5);
                REQUIRE(loaded_path_index.get_pangenome_pos("5", 11) == 12);
                REQUIRE(loaded_path_index.get_pangenome_pos("5-", 0) == 0);
                REQUIRE(loaded_path_index.get_pangenome_pos("5-", 1) == 1);
                REQUIRE(loaded_path_index.get_pangenome_pos("5-", 5) == 8);
                REQUIRE(loaded_path_index.get_pangenome_pos("5-", 12) == 6);
                REQUIRE(loaded_path_index.get_pangenome_pos("5-", 4) == 7);
                REQUIRE(loaded_path_index.get_pangenome_pos("5-", 11) == 5);
                REQUIRE(loaded_path_index.get_pangenome_pos("5-m", 0) == 0);
                REQUIRE(loaded_path_index.get_pangenome_pos("5-m", 1) == 1);
                REQUIRE(loaded_path_index.get_pangenome_pos("5-m", 5) == 8);
                REQUIRE(loaded_path_index.get_pangenome_pos("5-m", 12) == 5);
                REQUIRE(loaded_path_index.get_pangenome_pos("5-m", 4) == 7);
                REQUIRE(loaded_path_index.get_pangenome_pos("5-m", 11) == 6);

                /// The program will exit(1) here.
                // REQUIRE(loaded_path_index.get_pangenome_pos("5-", 24) == 0);
                // REQUIRE(loaded_path_index.get_pangenome_pos("4", 1) == 0);
                // REQUIRE(loaded_path_index.get_pangenome_pos("5", 24) == 0);
                // REQUIRE(loaded_path_index.get_pangenome_pos("4", 1) == 0);
            }
        }
    }
}
