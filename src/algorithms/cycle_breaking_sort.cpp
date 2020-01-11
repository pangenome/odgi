#include "cycle_breaking_sort.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

std::vector<handle_t> cycle_breaking_sort(const graph_t& graph) {
    std::vector<std::vector<uint64_t> > rank(graph.get_node_count());
    uint64_t i = 0;
    uint64_t j = 0;
    dfs(graph,
        [&](const handle_t& handle) { },
        [&](const handle_t& handle) {
            rank[number_bool_packing::unpack_number(handle)] = { j, i++, number_bool_packing::unpack_number(handle) };
        },
        [&](const handle_t& handle) { return false; },
        [&](void) { return false; },
        [&](const edge_t& e) { },
        [&](const edge_t& e) { ++j; },
        [&](const edge_t& e) { },
        [&](const edge_t& e) { },
        { },
        { });
    //std::cerr << "order size " << order.size() << " vs " <<  graph.node_size() << std::endl;
    std::sort(rank.begin(), rank.end());
    std::vector<handle_t> order; order.reserve(rank.size());
    for (auto& g : rank) {
        order.push_back(number_bool_packing::pack(g[2], false));
    }
    return order;
    //return order;
    /*
    graph_t tmp = graph;
    tmp.apply_ordering(order);
    return topological_order(&tmp);
    */
    /*
    // remove any edge whose from has a higher index than its to
    graph_t broken = graph;
    graph.for_each_edge([&](const edge_t& e) {
            if (graph.get_id(e.first) > graph.get_id(e.second)) {
                broken.destroy_edge(e);
                broken.create_edge(graph.get_is_reverse(e.second)?graph.flip(e.second):e.second,
                                   graph.get_is_reverse(e.first)?graph.flip(e.first):e.first);
            } else
            if (graph.get_is_reverse(e.first)
                       || graph.get_is_reverse(e.second)) {
                // force it into the graph as if it's not inverting
                broken.destroy_edge(e);
                broken.create_edge(graph.get_is_reverse(e.first)?graph.flip(e.first):e.first,
                                   graph.get_is_reverse(e.second)?graph.flip(e.second):e.second);
            }
            return true;
        });
    return topological_order(&broken);
    */
    /*
    // TODO... dagify, then obtain order, project back to previous graph
    graph_t tmp;
    ska::flat_hash_map<handlegraph::nid_t, handlegraph::nid_t> dag_to_base = dagify(&broken, &tmp, 0);
    std::vector<handlegraph::nid_t> pre_order;
    tmp.for_each_handle([&](const handle_t& handle) {
            pre_order.push_back(tmp.get_id(handle));
        });
    std::vector<handle_t> dag_order = lazy_topological_order(&tmp);
    ska::flat_hash_map<handlegraph::nid_t, std::vector<handlegraph::nid_t> > base_to_dag;
    //std::vector<handle_t> final_order;
    // take the order given by the first instance of the mapping to a particular handle in the base graph
    for (auto& handle : dag_order) {
        //tmp.for_each_handle([&](const handle_t& handle) {
        handlegraph::nid_t dag_id = tmp.get_id(handle);
        auto& base_id = dag_to_base[dag_id];
        //if (base_to_dag.find(base_id) == base_to_dag.end()) {
        base_to_dag[base_id].push_back(dag_id);
        //final_order.push_back(graph.get_handle(base_id));
    }
        //});
    std::vector<std::pair<handlegraph::nid_t, handlegraph::nid_t> > x;
    graph.for_each_handle([&](const handle_t& handle) {
            handlegraph::nid_t base_id = graph.get_id(handle);
            //auto& dag_id = base_to_dag[base_id].front();
            auto& dag_v = base_to_dag[base_id];
            //std::sort(dag_v.begin(), dag_v.end());
            x.push_back(make_pair(dag_v.back(), base_id));
            //final_order.push_back(graph.get_handle(base_id));
        });
    std::sort(x.begin(), x.end());
    std::vector<handle_t> final_order;
    for (auto p : x) {
        final_order.push_back(graph.get_handle(p.second));
    }
    return final_order;
    */
}

}

}
