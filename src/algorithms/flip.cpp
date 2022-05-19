#include "flip.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

void flip_paths(graph_t& graph,
                const std::vector<path_handle_t>& no_flips) {
    // for each path, find its average orientation
    ska::flat_hash_set<path_handle_t> no_flip;
    for (auto& p : no_flips) { no_flip.insert(p); }
    std::vector<path_handle_t> paths;
    graph.for_each_path_handle([&](const path_handle_t& p) { paths.push_back(p); });
    ska::flat_hash_map<path_handle_t, path_handle_t> to_flip;
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
                auto flipped = graph.create_path_handle(name);
#pragma omp critical (to_flip)
                to_flip[path] = flipped;
            }
        }
    }
    // for each, add a new path named "_inv" and write it in the reverse orientation
    std::vector<path_handle_t> to_flip_v;
    for (auto& p : to_flip) { to_flip_v.push_back(p.first); }
#pragma omp parallel for
    for (auto& path : to_flip_v) {
        auto& flipped = to_flip[path];
        std::vector<handle_t> v;
        graph.for_each_step_in_path(
            path,
            [&v,&graph](const step_handle_t& s) {
                auto h = graph.flip(graph.get_handle_of_step(s));
                v.push_back(h);
            });
        std::reverse(v.begin(), v.end());
        for (auto& h : v) {
            graph.append_step(flipped, h);
        }
    }
    // apply original order and record which paths to drop
    std::vector<path_handle_t> order;
    std::vector<path_handle_t> to_drop;
    for (auto& path : paths) {
        auto f = to_flip.find(path);
        if (f != to_flip.end()) {
            order.push_back(path);
            to_drop.push_back(as_path_handle(order.size()));
            order.push_back(f->second);
        } else {
            order.push_back(path);
        }
    }
    // set the order
    graph.apply_path_ordering(order);
    // delete the old paths, then optimize the graph and return it (by reference)
#pragma omp parallel for
    for (auto& path : to_drop) {
        graph.destroy_path(path);
    }
    // clean up deleted path steps
    graph.optimize();
}

}

}
