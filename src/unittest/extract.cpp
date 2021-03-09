#include "catch.hpp"

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include "odgi.hpp"


#include <iostream>
#include <random>

#include "src/algorithms/subgraph/extract.hpp"
#include "src/algorithms/subgraph/region.hpp"

namespace odgi {

    namespace unittest {

    using namespace std;
    using namespace handlegraph;

        TEST_CASE("Extracting a path range", "[extracting]") {
            graph_t graph;
            handle_t n1 = graph.create_handle("CAAATAAG");
            handle_t n2 = graph.create_handle("A");
            handle_t n3 = graph.create_handle("G");
            handle_t n4 = graph.create_handle("T");
            handle_t n5 = graph.create_handle("C");
            handle_t n6 = graph.create_handle("TTG");
            graph.create_edge(n1, n2);
            graph.create_edge(n1, n3);
            graph.create_edge(n2, n3);
            graph.create_edge(n2, n4);
            graph.create_edge(n3, n4);
            graph.create_edge(n3, n5);
            graph.create_edge(n4, n6);
            graph.create_edge(n5, n6);

            auto path_x = graph.create_path_handle("x");
            graph.append_step(path_x, n1);
            graph.append_step(path_x, n3);
            graph.append_step(path_x, n5);
            graph.append_step(path_x, n6);

            auto path_y = graph.create_path_handle("y");
            graph.append_step(path_y, n1);
            graph.append_step(path_y, n2);
            graph.append_step(path_y, n5);
            graph.append_step(path_y, n4);

            auto path_z = graph.create_path_handle("z");
            graph.append_step(path_z, n1);
            graph.append_step(path_z, n3);
            graph.append_step(path_z, n3);
            graph.append_step(path_z, n6);

            Region region {"x", 8, 10};
            uint64_t context_size = 0;

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

            graph_t subgraph;
            algorithms::extract_path_range(graph, path_x, region.start, region.end, subgraph);
            algorithms::expand_subgraph_by_steps(graph, subgraph, context_size, false);
            algorithms::add_connecting_edges_to_subgraph(graph, subgraph);
            algorithms::add_subpaths_to_subgraph(graph, subgraph, false);

            subgraph.for_each_path_handle(test_path);

            SECTION("The extracted subgraph is as expected") {
                REQUIRE(subgraph.get_path_count() == 3);

                REQUIRE(subgraph.get_node_count() == 3);
                REQUIRE(subgraph.get_sequence(subgraph.get_handle(3)) == "G");
                REQUIRE(subgraph.get_sequence(subgraph.get_handle(5)) == "C");
                REQUIRE(subgraph.get_sequence(subgraph.get_handle(6)) == "TTG");

                auto new_path_x = subgraph.get_path_handle("x");
                REQUIRE(subgraph.get_handle_of_step(subgraph.path_begin(new_path_x)) == graph.get_handle(3));
                std::string x;
                subgraph.for_each_step_in_path(new_path_x, [&](const step_handle_t& step) {
                    x.append(graph.get_sequence(subgraph.get_handle_of_step(step)));
                });
                REQUIRE(x == "GCTTG");

                auto new_path_y = subgraph.get_path_handle("y");
                REQUIRE(subgraph.get_handle_of_step(subgraph.path_begin(new_path_y)) == graph.get_handle(5));
                std::string y;
                subgraph.for_each_step_in_path(new_path_y, [&](const step_handle_t& step) {
                    y.append(graph.get_sequence(subgraph.get_handle_of_step(step)));
                });
                REQUIRE(y == "C");

                auto new_path_z = subgraph.get_path_handle("z");
                REQUIRE(subgraph.get_handle_of_step(subgraph.path_begin(new_path_z)) == graph.get_handle(3));
                std::string z;
                subgraph.for_each_step_in_path(new_path_z, [&](const step_handle_t& step) {
                    z.append(graph.get_sequence(subgraph.get_handle_of_step(step)));
                });
                REQUIRE(z == "GGTTG");
            }

        }

    }
}
