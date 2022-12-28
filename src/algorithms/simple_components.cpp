#include "simple_components.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

// the set of components that could be merged into single nodes without changing the path space of the graph
std::vector<std::vector<handle_t>> simple_components(
    // todo unused return_all_handles
    const PathHandleGraph &graph, const uint64_t& min_size, const bool& return_all_handles, const uint64_t& nthreads) {

    std::vector<uint64_t> data; data.reserve(graph.get_node_count()*2);
    graph.for_each_handle(
        [&](const handle_t& handle) {
            data.push_back(as_integer(handle));
            data.push_back(as_integer(graph.flip(handle)));
        });
    boophf_uint64_t bphf(data.size(),data,nthreads,2.0,false,false); // mapping structure for original node ids

    // todo check if we should or shouldn't use the gcc atomic primitives
    std::vector<std::atomic<DisjointSets::Aint>> simple_data(data.size()+1); // maps into this set of disjoint sets
    auto simple_dset = DisjointSets(simple_data.data(), simple_data.size());

    bool in_parallel = nthreads > 1;

    graph.for_each_handle(
        [&](const handle_t& h) {
            //uint64_t h_i = bphf.lookup(as_integer(h));
            //uint64_t h_j = bphf.lookup(as_integer(graph.flip(h)));
            if (graph.get_degree(h, true) == 1) {
                // go backward
                graph.follow_edges(
                    h, true,
                    [&](const handle_t& prev) {
                        if (graph.get_id(h) != graph.get_id(prev)
                            && graph.get_degree(prev, false) == 1
                            && nodes_are_perfect_path_neighbors(graph, prev, h)) {
                            uint64_t from_i = bphf.lookup(as_integer(prev));
                            uint64_t to_i = bphf.lookup(as_integer(h));
                            simple_dset.unite(from_i, to_i);
                        }
                    });
            }
            if (graph.get_degree(h, false) == 1) {
                // go forward
                graph.follow_edges(
                    h, false,
                    [&](const handle_t& next) {
                        if (graph.get_id(h) != graph.get_id(next)
                            && graph.get_degree(next, true) == 1
                            && nodes_are_perfect_path_neighbors(graph, h, next)) {
                            uint64_t from_i = bphf.lookup(as_integer(h));
                            uint64_t to_i = bphf.lookup(as_integer(next));
                            simple_dset.unite(from_i, to_i);
                        }
                    });
            }
        }, in_parallel);

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
    auto pred_in_comp = [&](const handle_t& h,
                            const ska::flat_hash_set<uint64_t>& comp) {
                            bool in = true;
                            graph.follow_edges(h, true, [&](const handle_t& prev) {
                                                            in &= comp.count(graph.get_id(prev));
                                                            //in &= graph.get_id(h) != graph.get_id(prev);
                                                        });
                            return in;
                        };

    for (auto& c : simple_components) {
        auto& comp = c.second;
        assert(comp.size());
        /*
        std::cerr << "on component " << i++ << std::endl;
        for (auto& h : comp) {
            std::cerr << graph.get_id(h) << (graph.get_is_reverse(h) ? "-" : "+") << ",";
        }
        std::cerr << std::endl;
        */

        if (comp.size() >= min_size) {
            // start somewhere
            ska::flat_hash_set<uint64_t> comp_set;
            for (auto& h : comp) comp_set.insert(graph.get_id(h));
            auto h_itr = comp.begin();
            // try to find the start of the component
            while (h_itr != comp.end()
                   && graph.get_degree(*h_itr, true) == 1
                   && pred_in_comp(*h_itr, comp_set)) {
                ++h_itr;
            }
            handle_t h;
            if (h_itr == comp.end()) {
                // could be looping, take the start as the first node in the component
                h = comp.front();
            } else {
                h = *h_itr;
            }
            // To avoid pushing the same handle (h in sorted_comp) until you have comp.size() elements
            if (graph.get_degree(h, false) > 0) {
                handle_components.emplace_back();
                // walk from our start node through the component
                auto& sorted_comp = handle_components.back();
                bool fail = false;
                do {
                    sorted_comp.push_back(h);
                    graph.follow_edges(h, false, [&](const handle_t& next) {
                                                    if (comp_set.count(graph.get_id(next))) {
                                                        h = next;
                                                    } else {
                                                        // previous ordering assumptions are violated
                                                        // this can be caused by self-looping head or tail nodes
                                                        // in the component
                                                        fail = true;
                                                    }
                                                });
                } while (sorted_comp.size() < comp.size() && !fail);
                if (sorted_comp.size() < min_size
                    || sorted_comp.size() < comp.size()) {
                    handle_components.pop_back();
                } else {
                    assert(sorted_comp.front() != sorted_comp.back());
                }
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
