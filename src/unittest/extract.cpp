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
            graph.create_handle("CAAA");
            graph.create_handle("AT");
            graph.create_handle("GCC");
            graph.create_handle("TA");
            graph.create_handle("CTT");
            graph.create_handle("TTGA");

            auto path_x = graph.create_path_handle("x");
            graph.append_step(path_x, graph.get_handle(1, true));
            graph.append_step(path_x, graph.get_handle(3, true));
            graph.append_step(path_x, graph.get_handle(2, true));
            graph.append_step(path_x, graph.get_handle(5, true));
            graph.append_step(path_x, graph.get_handle(6, true));

            auto path_y = graph.create_path_handle("y");
            graph.append_step(path_y, graph.get_handle(1, false));
            graph.append_step(path_y, graph.get_handle(2, true));
            graph.append_step(path_y, graph.get_handle(4, false));
            graph.append_step(path_y, graph.get_handle(5, false));
            graph.append_step(path_y, graph.get_handle(4, false));

            auto path_z = graph.create_path_handle("z");
            graph.append_step(path_z, graph.get_handle(1, false));
            graph.append_step(path_z, graph.get_handle(3, false));
            graph.append_step(path_z, graph.get_handle(3, true));
            graph.append_step(path_z, graph.get_handle(6, false));

            Region region {"x", 4, 8};
            uint64_t context_size = 0;

            graph_t subgraph;
            algorithms::extract_path_range(graph, path_x, region.start, region.end, subgraph);
            algorithms::expand_subgraph_by_steps(graph, subgraph, context_size, false);
            algorithms::add_connecting_edges_to_subgraph(graph, subgraph);

            std::vector<path_handle_t> paths;
            paths.reserve(graph.get_path_count());
            graph.for_each_path_handle([&](const path_handle_t path) {
                paths.push_back(path);
            });
            algorithms::add_subpaths_to_subgraph(graph, paths, subgraph, 2);

            auto test_path = [&](const path_handle_t& p) {
                auto& path_meta = subgraph.path_metadata(p);
                uint64_t i = 0;
                subgraph.for_each_step_in_path(p, [&](const step_handle_t& step) {
                    handle_t h = subgraph.get_handle_of_step(step);
                    ++i;
                });
                REQUIRE(i == path_meta.length);
            };

            subgraph.for_each_path_handle(test_path);

            SECTION("The extracted subgraph is as expected") {
                REQUIRE(subgraph.get_path_count() == 3);

                REQUIRE(subgraph.get_node_count() == 2);
                REQUIRE(subgraph.get_sequence(subgraph.get_handle(2)) == "AT");
                REQUIRE(subgraph.get_sequence(subgraph.get_handle(3)) == "GCC");

                REQUIRE(subgraph.has_path("x:4-9"));
                auto new_path_x = subgraph.get_path_handle("x:4-9");
                REQUIRE(subgraph.get_handle_of_step(subgraph.path_begin(new_path_x)) == graph.get_handle(3, true));
                std::string x;
                subgraph.for_each_step_in_path(new_path_x, [&](const step_handle_t& step) {
                    x.append(subgraph.get_sequence(subgraph.get_handle_of_step(step)));
                });
                REQUIRE(x == "GGCAT");

                REQUIRE(subgraph.has_path("y:4-6"));
                auto new_path_y = subgraph.get_path_handle("y:4-6");
                REQUIRE(subgraph.get_handle_of_step(subgraph.path_begin(new_path_y)) == graph.get_handle(2, true));
                std::string y;
                subgraph.for_each_step_in_path(new_path_y, [&](const step_handle_t& step) {
                    y.append(subgraph.get_sequence(subgraph.get_handle_of_step(step)));
                });
                REQUIRE(y == "AT");

                REQUIRE(subgraph.has_path("z:4-10"));
                auto new_path_z = subgraph.get_path_handle("z:4-10");
                REQUIRE(subgraph.get_handle_of_step(subgraph.path_begin(new_path_z)) == graph.get_handle(3));
                std::string z;
                subgraph.for_each_step_in_path(new_path_z, [&](const step_handle_t& step) {
                    z.append(subgraph.get_sequence(subgraph.get_handle_of_step(step)));
                });
                REQUIRE(z == "GCCGGC");
            }

        }

    }
}
