#include "flip.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

void flip_paths(graph_t& graph,
                graph_t& into,
                const std::vector<path_handle_t>& no_flips,
                const std::vector<path_handle_t>& ref_flips) {
    // Copy the nodes
    graph.for_each_handle([&](const handle_t& h) {
        into.create_handle(graph.get_sequence(h), graph.get_id(h));
    });
    // Copy the edges
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

    // Paths to not flip
    ska::flat_hash_set<path_handle_t> no_flip;
    for (auto& p : no_flips) { no_flip.insert(p); }

    // Paths to use as reference for flipping
    ska::flat_hash_set<path_handle_t> ref_flip;
    for (auto& p : ref_flips) { ref_flip.insert(p); }

    // Prepare all path handles in a vector (for parallel processing)
    std::vector<path_handle_t> paths;
    // for each path, find its average orientation
    graph.for_each_path_handle([&](const path_handle_t& p) { paths.push_back(p); });

    // Check whether reference paths are (in general) in fwd or rev orientation
    uint64_t ref_rev = 0;
    uint64_t ref_fwd = 0;
    // OpenMP's reduction: each thread maintains its own private copy of these variables
    // during the parallel region, and OpenMP combines them at the end using addition.
#pragma omp parallel for reduction(+:ref_rev,ref_fwd)
    for (auto& path : paths) {
        if (ref_flip.count(path)) {
            graph.for_each_step_in_path(
                path,
                [&ref_rev,&ref_fwd,&graph](const step_handle_t& s) {
                    auto h = graph.get_handle_of_step(s);
                    auto len = graph.get_length(h);
                    if (graph.get_is_reverse(h)) {
                        ref_rev += len;
                    } else {
                        ref_fwd += len;
                    }
                });
        }
    }

#pragma omp parallel for
    for (auto& path : paths) {
        // ref_paths must not be flipped either
        if (!no_flip.count(path) && !ref_flip.count(path)) {
            // Check path orientation with respect to the graph
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

            // Check if the path should be flipped
            bool flip_path = false;
            if (ref_flip.size() > 0) {
                // if ref paths are reversed, reversed paths should be
                // not flipped to stay consistent with the reference
                flip_path = ref_rev > ref_fwd ? rev < fwd : rev > fwd;
            } else {
                // those that tend to be reversed more than forward should be flipped
                flip_path = rev > fwd;
            }
            if (flip_path) {
                auto name = graph.get_path_name(path) + "_inv";
                auto flipped = into.create_path_handle(name);
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
