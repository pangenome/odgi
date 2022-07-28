#include "flip.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

void flip_paths(graph_t& graph,
                graph_t& into,
                const std::vector<path_handle_t>& no_flips) {
    // for each path, find its average orientation
    graph.for_each_handle([&](const handle_t& h) {
        into.create_handle(graph.get_sequence(h), graph.get_id(h));
    });
    graph.for_each_handle([&](const handle_t& h) {
        graph.follow_edges(h, false, [&](const handle_t& t) {
            into.create_edge(into.get_handle(graph.get_id(h), graph.get_is_reverse(h)),
                             into.get_handle(graph.get_id(t), graph.get_is_reverse(t)));
        });
        auto r = graph.flip(h);
        graph.follow_edges(r, false, [&](const handle_t& t) {
            into.create_edge(into.get_handle(graph.get_id(r), graph.get_is_reverse(r)),
                             into.get_handle(graph.get_id(t), graph.get_is_reverse(t)));
        });
    });
    ska::flat_hash_set<path_handle_t> no_flip;
    for (auto& p : no_flips) { no_flip.insert(p); }
    std::vector<path_handle_t> paths;
    std::vector<path_handle_t> into_paths;
    graph.for_each_path_handle([&](const path_handle_t& p) { paths.push_back(p); });
#pragma omp parallel for
    for (auto& path : paths) {
        if (!no_flip.count(path)) {
            uint64_t rev = 0;
            uint64_t fwd = 0;
            graph.for_each_step_in_path(
                path,
                [&rev,&fwd,&graph](const step_handle_t& s) {
                    auto h = graph.get_handle_of_step(s);
                    auto len = graph.get_length(h);
                    if (graph.get_is_reverse(h)) {
                        rev += len;
                    } else {
                        fwd += len;
                    }
                });
            // those that tend to be reversed more than forward should be flipped
            if (rev > fwd) {
                auto name = graph.get_path_name(path) + "_inv";
                auto flipped = into.create_path_handle(name);
                into_paths.push_back(flipped);
                std::vector<handle_t> v;
                graph.for_each_step_in_path(
                    path,
                    [&into,&v,&graph](const step_handle_t& s) {
                        auto h = graph.flip(graph.get_handle_of_step(s));
                        auto q = into.get_handle(graph.get_id(h), graph.get_is_reverse(h));
                        v.push_back(q);
                    });
                std::reverse(v.begin(), v.end());
                for (auto& q : v) {
                    into.append_step(flipped, q);
                }
            } else {
                auto fwd = into.create_path_handle(graph.get_path_name(path));
                into_paths.push_back(fwd);
                graph.for_each_step_in_path(
                    path,
                    [&fwd,&into,&graph](const step_handle_t& s) {
                        auto h = graph.get_handle_of_step(s);
                        auto q = into.get_handle(graph.get_id(h), graph.get_is_reverse(h));
                        into.append_step(fwd, q);
                    });
            }
        } else {
            // add the path as-is
            auto fwd = into.create_path_handle(graph.get_path_name(path));
            into_paths.push_back(fwd);
            graph.for_each_step_in_path(
                path,
                [&fwd,&into,&graph](const step_handle_t& s) {
                    auto h = graph.get_handle_of_step(s);
                    auto q = into.get_handle(graph.get_id(h), graph.get_is_reverse(h));
                    into.append_step(fwd, q);
                });
        }
    }
    into.optimize();
}

}

}
