#include "flip.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

void flip_paths(graph_t& graph,
                graph_t& into,
                const std::vector<path_handle_t>& no_flips) {
    graph.for_each_handle([&](const handle_t& h) {
        into.create_handle(graph.get_sequence(h), graph.get_id(h));
    });
    graph.for_each_handle([&](const handle_t& h) {
        graph.follow_edges(h, false, [&](const handle_t& next) {
            into.create_edge(into.get_handle(graph.get_id(h), graph.get_is_reverse(h)),
                             into.get_handle(graph.get_id(next), graph.get_is_reverse(next)));
        });
        graph.follow_edges(h, true, [&](const handle_t& prev) {
            into.create_edge(into.get_handle(graph.get_id(prev), graph.get_is_reverse(prev)),
                             into.get_handle(graph.get_id(h), graph.get_is_reverse(h)));
        });
    });
    ska::flat_hash_set<path_handle_t> no_flip;
    for (auto& p : no_flips) { no_flip.insert(p); }
    std::vector<path_handle_t> paths;
    std::vector<path_handle_t> into_paths;
    // for each path, find its average orientation
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

    ska::flat_hash_set<std::pair<handle_t, handle_t>> edges_to_create;

#pragma omp parallel for
    for (auto& path : paths) {
        // New edges can be due only when paths are flipped
        if (!no_flip.count(path)) {
            handle_t last;
            const step_handle_t begin_step = into.path_begin(path);
            into.for_each_step_in_path(path, [&](const step_handle_t &step) {
                handle_t h = into.get_handle_of_step(step);
                if (step != begin_step && !into.has_edge(last, h)) {
    #pragma omp critical (edges_to_create)
                    edges_to_create.insert({last, h});
                }
                last = h;
            });
        }
    }

    // add missing edges
    for (auto edge: edges_to_create) {
        into.create_edge(edge.first, edge.second);
    }

    into.optimize();
}

}

}
