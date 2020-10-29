#include "simple_components.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

// the set of components that could be merged into single nodes without changing the path space of the graph
std::vector<std::vector<handle_t>> simple_components(
    const PathHandleGraph &graph, const uint64_t& min_size, const bool& return_all_handles, const uint64_t& nthreads) {

    std::vector<uint64_t> data; data.reserve(graph.get_node_count()*2);
    graph.for_each_handle(
        [&](const handle_t& handle) {
            data.push_back(as_integer(handle));
            data.push_back(as_integer(graph.flip(handle)));
        });
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

    bool in_parallel = nthreads > 1;
    uint64_t curr_threads = get_thread_count();
    if (in_parallel) omp_set_num_threads(nthreads);
    graph.for_each_edge(
        [&](const edge_t& edge) {
            const handle_t& from = edge.first;
            const handle_t& to = edge.second;
            // unite the handles with themselves
            self_unite(from);
            self_unite(to);
            if (graph.get_id(from) != graph.get_id(to)) {
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
            }
        }, in_parallel);
    if (in_parallel) omp_set_num_threads(curr_threads);

    ska::flat_hash_map<uint64_t, std::vector<handle_t>> simple_components;
    graph.for_each_handle(
        [&](const handle_t& handle) {
            uint64_t a_id = simple_dset.find(bphf.lookup(as_integer(handle)));
            // take the component of our forward orientation
            simple_components[a_id].push_back(handle);
        });

    // now we combine and order the nodes in each dset
    std::vector<std::vector<handle_t>> handle_components;
    //std::cerr << "processing components" << std::endl;
    for (auto& c : simple_components) {
        //std::cerr << "on component " << i++ << std::endl;
        auto& comp = c.second;
        std::sort(
            comp.begin(), comp.end(),
            [&](const handle_t& a, const handle_t& b) {
                return as_integer(a) > as_integer(b);
            });
        assert(comp.size());
        if (comp.size() >= min_size) {
            // start somewhere
            ska::flat_hash_set<uint64_t> comp_set;
            for (auto& h : comp) comp_set.insert(as_integer(h));
            handle_t h = comp.front();
            handle_t base = h; // so we don't loop endlessly in self-linking components
            //std::cerr << "base is " << graph.get_id(base) << std::endl;
            bool has_prev = false;
            do {
                has_prev = graph.get_degree(h, true) == 1;
                handle_t prev = h;
                if (has_prev) {
                    graph.follow_edges(h, true, [&](const handle_t& p) { prev = p; });
                }
                if (h != prev
                    && prev != base
                    && comp_set.count(as_integer(prev))) {
                    //std::cerr << "continuing to prev = " << graph.get_id(prev) << std::endl;
                    h = prev;
                } else {
                    //std::cerr << "breaking at prev = " << graph.get_id(prev) << std::endl;
                    has_prev = false;
                }
            } while (has_prev);
            base = h; // reset base
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
                if (h != next
                    && next != base
                    && comp_set.count(as_integer(next))) {
                    //std::cerr << "continuing to next = " << graph.get_id(next) << std::endl;
                    h = next;
                } else {
                    //std::cerr << "breaking at next = " << graph.get_id(next) << std::endl;
                    has_next = false;
                }
            } while (has_next);
            assert(sorted_comp.size() > 0);
            if (sorted_comp.size() < min_size) {
                handle_components.pop_back();
            }
        }
    }

    if (min_size == 1) {
        uint64_t seen_nodes = 0;
        for (auto& comp : simple_components) {
            seen_nodes += comp.second.size();
        }
        assert(seen_nodes == graph.get_node_count());
    }

    return handle_components;
}
}

}
