#include "catch.hpp"

#include "algorithms/cut_tips.hpp"
#include "odgi.hpp"

#include <vector>

namespace odgi {
namespace unittest {

using namespace handlegraph;

TEST_CASE("Copying a graph preserves path traversal", "[paths][copy][regression]") {
    graph_t source;
    auto removed_path = source.create_path_handle("removed");
    source.destroy_path(removed_path);

    std::vector<handle_t> source_handles;
    for (const auto& sequence : {"A", "C", "G", "T", "N"}) {
        source_handles.push_back(source.create_handle(sequence));
    }

    auto source_path = source.create_path_handle("path");
    for (const auto& handle : source_handles) {
        source.append_step(source_path, handle);
    }
    REQUIRE(source.get_id(source.get_handle_of_step(source.path_begin(source_path))) == source.get_id(source_handles.front()));

    graph_t copy;
    copy.copy(source);

    auto copied_path = copy.get_path_handle("path");
    auto step = copy.path_begin(copied_path);
    for (const auto& handle : source_handles) {
        REQUIRE(copy.get_id(copy.get_handle_of_step(step)) == source.get_id(handle));
        step = copy.get_next_step(step);
    }
    REQUIRE(step == copy.path_end(copied_path));
}

TEST_CASE("Cutting boundary tips preserves the remaining path", "[paths][cut_tips][regression]") {
    graph_t graph;
    auto left = graph.create_handle("A");
    auto middle = graph.create_handle("C");
    auto right = graph.create_handle("G");
    graph.create_edge(left, middle);
    graph.create_edge(middle, right);

    auto path = graph.create_path_handle("path");
    graph.append_step(path, left);
    graph.append_step(path, middle);
    graph.append_step(path, right);

    algorithms::cut_tips(graph, 0);

    REQUIRE(graph.get_step_count(path) == 1);
    auto remaining = graph.path_begin(path);
    REQUIRE(graph.get_id(graph.get_handle_of_step(remaining)) == graph.get_id(middle));
    REQUIRE_FALSE(graph.has_previous_step(remaining));
    REQUIRE_FALSE(graph.has_next_step(remaining));

    graph.optimize();
    std::vector<nid_t> path_ids;
    graph.for_each_step_in_path(path, [&](const step_handle_t& step) {
        path_ids.push_back(graph.get_id(graph.get_handle_of_step(step)));
    });
    REQUIRE(path_ids == std::vector<nid_t>{graph.get_id(graph.get_handle(1))});
}

}
}
