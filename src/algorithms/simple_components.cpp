#include "simple_components.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

// the set of components that could be merged into single nodes without changing the path space of the graph
std::vector<std::vector<handle_t>> simple_components(
    const PathHandleGraph &graph, const uint64_t& min_size, const bool return_all_handles) {

    std::vector<uint64_t> data; data.reserve(graph.get_node_count()*2);
    graph.for_each_handle(
        [&](const handle_t& handle) {
            data.push_back(as_integer(handle));
            data.push_back(as_integer(graph.flip(handle)));
        });
    uint64_t nthreads = get_thread_count();
    boophf_uint64_t bphf(data.size(),data,nthreads,2.0,false,false);

    // todo check if we should or shouldn't use the gcc atomic primitives
    std::vector<DisjointSets::Aint> simple_data(data.size()+1);
    auto simple_dset = DisjointSets(simple_data.data(), simple_data.size());

    auto self_unite =
        [&](const handle_t& h) {
            uint64_t h_i = bphf.lookup(as_integer(h));
            uint64_t h_j = bphf.lookup(as_integer(graph.flip(h)));
            simple_dset.unite(h_i, h_j);
        };

    graph.for_each_edge(
        [&](const edge_t& edge) {
            const handle_t& from = edge.first;
            const handle_t& to = edge.second;
            // unite the handles with themselves
            self_unite(from);
            self_unite(to);
            // check if they can be linked across the edge
            if (graph.get_degree(from, false) == 1
                && graph.get_degree(to, true) == 1
                && nodes_are_perfect_path_neighbors(graph, from, to)) {
                //std::cerr << "would merge " << graph.get_id(from)
                //          << " with " << graph.get_id(to) << std::endl;
                uint64_t from_i = bphf.lookup(as_integer(from));
                uint64_t to_i = bphf.lookup(as_integer(to));
                //std::cerr << "indexes " << from_i << " and " << to_i << std::endl;
                simple_dset.unite(from_i, to_i);
            }
        },
        true); // parallel

    ska::flat_hash_map<uint64_t, std::vector<handle_t>> simple_components;
    graph.for_each_handle(
        [&](const handle_t& handle) {
            uint64_t a_id = simple_dset.find(bphf.lookup(as_integer(handle)));
            // take the component of our forward orientation
            simple_components[a_id].push_back(handle);
        });

    // now we combine and order the nodes in each dset
    std::vector<std::vector<handle_t>> handle_components;
    for (auto& c : simple_components) {
        auto& comp = c.second;
        if (comp.size() >= min_size) {
            // start somewhere
            ska::flat_hash_set<uint64_t> comp_set;
            for (auto& h : comp) comp_set.insert(as_integer(h));
            handle_t h = comp.front();
            bool has_prev = false;
            do {
                has_prev = graph.get_degree(h, true) == 1;
                handle_t prev = h;
                if (has_prev) {
                    graph.follow_edges(h, true, [&](const handle_t& p) { prev = p; });
                }
                if (comp_set.count(as_integer(prev))) {
                    h = prev;
                } else {
                    has_prev = false;
                }
            } while (has_prev);

            //std::cerr << "Front is " << graph.get_id(h) << std::endl;
            handle_components.emplace_back();
            auto& sorted_comp = handle_components.back();
            bool has_next = false;
            do {
                sorted_comp.push_back(h);
                has_next = graph.get_degree(h, false) == 1;
                handle_t next = h;
                if (has_next) {
                    graph.follow_edges(h, false, [&](const handle_t& n) { next = n; });
                }
                if (comp_set.count(as_integer(next))) {
                    h = next;
                } else {
                    has_next = false;
                }
            } while (has_next);
        }
    }

    return handle_components;
}
}

}
