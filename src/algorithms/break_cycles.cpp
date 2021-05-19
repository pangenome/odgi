#include "break_cycles.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

std::vector<edge_t> edges_inducing_cycles(
    const HandleGraph& graph,
    const uint64_t& max_cycle_size,
    const uint64_t& max_search_bp
    ) {

    uint64_t min_handle_rank = 0;
    uint64_t max_handle_rank = 0;
    graph.for_each_handle([&](const handle_t& found) {
        uint64_t handle_rank = number_bool_packing::unpack_number(found);
        min_handle_rank = std::min(min_handle_rank, handle_rank);
        max_handle_rank = std::max(max_handle_rank, handle_rank);
    });

    handle_t min_handle = number_bool_packing::pack(min_handle_rank, false);
    handle_t root_handle = min_handle;
    //std::vector<handle_t> seeds = { min_handle };

    // edges that we'll find greedily and which we should remove
    ska::flat_hash_set<edge_t> edges_to_remove;

    graph.for_each_handle(
        [&](const handle_t& handle) {
            //std::cerr << "on handle " << graph.get_id(root_handle) << std::endl;
            uint64_t seen_bp = 0;
            uint64_t max_depth = 0;
            uint64_t last_min_length_bp = 0; //std::numeric_limits<uint64_t>::max();
            uint64_t curr_min_length_bp = std::numeric_limits<uint64_t>::max();
            handle_t handle_fwd = handle;
            handle_t handle_rev = graph.flip(handle);
            for (const handle_t& root_handle : { handle_fwd, handle_rev }) {
                bfs(graph,
                    [&graph,&edges_to_remove,&seen_bp,&curr_min_length_bp,&last_min_length_bp,&max_depth]
                    (const handle_t& h, const uint64_t& r, const uint64_t& l, const uint64_t& d) {
                        /*
                          std::cerr << "stepping onto " << graph.get_id(h) << (graph.get_is_reverse(h)?"-":"+")
                          << " (" << r << " " << l << " " << d << ")"
                          << std::endl;
                        */
                        // track the minimum length at the last step
                        if (d > max_depth) {
                            max_depth = d;
                            last_min_length_bp = curr_min_length_bp;
                            curr_min_length_bp = l;
                        } else {
                            curr_min_length_bp = std::min(l, curr_min_length_bp);
                        }
                        seen_bp += graph.get_length(h);
                    },
                    [](const handle_t& h) { return false; },
                    [&graph,&root_handle,&edges_to_remove](const handle_t& p, const handle_t& h) {
                        edge_t e = std::make_pair(p, h);
                        if (h == root_handle) {
                            edges_to_remove.insert(e);
                            return true;
                        } else if (edges_to_remove.count(e)
                                   || edges_to_remove.count(std::make_pair(graph.flip(h), graph.flip(p)))) {
                            return true;
                        } else {
                            return false;
                        }
                    },
                    [&last_min_length_bp,&max_cycle_size,&seen_bp,&max_search_bp]() {
                        //std::cerr << "min_length_bp " << last_min_length_bp << " seen_bp " << seen_bp << std::endl;
                        return last_min_length_bp > max_cycle_size || seen_bp > max_search_bp;
                    },
                    { root_handle },
                    { },
                    false); // don't use bidirectional search
            }
        });

    std::vector<edge_t> edges;
    for (auto& e : edges_to_remove) {
        edges.push_back(e);
    }

    return edges;
}

std::vector<edge_t> edges_inducing_cycles_iter(
    const HandleGraph& graph,
    const uint64_t& max_cycle_size,
    const uint64_t& max_search_bp,
    const uint64_t& iter_max) {
    std::vector<edge_t> cycle_edges;
    for (uint64_t i = 0; i < iter_max; ++i) {
        std::vector<edge_t> new_cycle_edges
            = algorithms::edges_inducing_cycles(graph,
                                                max_cycle_size,
                                                max_search_bp);
        if (new_cycle_edges.empty()) break;
        cycle_edges.reserve(cycle_edges.size() + new_cycle_edges.size());
        cycle_edges.insert(cycle_edges.end(), new_cycle_edges.begin(), new_cycle_edges.end());
    }
    return cycle_edges;
}

uint64_t break_cycles(
    DeletableHandleGraph& graph,
    const uint64_t& max_cycle_size,
    const uint64_t& max_search_bp,
    const uint64_t& iter_max) {

    uint64_t removed_edges = 0;
    for (uint64_t i = 0; i < iter_max; ++i) {
        std::vector<edge_t> edges_to_remove
            = algorithms::edges_inducing_cycles(graph,
                                                max_cycle_size,
                                                max_search_bp);
        if (edges_to_remove.empty()) {
            break;
        }
        for (auto& edge : edges_to_remove) {
            graph.destroy_edge(edge);
        }
        removed_edges += edges_to_remove.size();
    }

    return removed_edges;
}

}

}
