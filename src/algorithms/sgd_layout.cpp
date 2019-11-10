#include "sgd_layout.hpp"
#include "sgd2.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

std::vector<double> sgd_layout(const HandleGraph& graph, uint64_t pivots, uint64_t t_max, double eps, double x_padding) {
    std::vector<double> layout(graph.get_node_count()*2);
    double max_x = 0;
    auto weak_components = algorithms::weakly_connected_components(&graph);
    for (auto& component : weak_components) {
        std::vector<handlegraph::nid_t> component_ids;
        ska::flat_hash_map<handlegraph::nid_t, uint64_t> local_id;
        for (auto& id : component) {
            local_id[id] = component_ids.size();
            component_ids.push_back(id);
        }
        // convert to input format for SGD
        std::vector<uint64_t> I, J;
        graph.for_each_edge([&](const edge_t& e) {
                if (component.count(graph.get_id(e.first))
                    && component.count(graph.get_id(e.first))) {
                    I.push_back(local_id[graph.get_id(e.first)]);
                    J.push_back(local_id[graph.get_id(e.second)]);
                }
            });
        uint64_t n = component.size();
        std::vector<double> X(2*n);
        std::random_device dev;
        // todo, seed with graph topology/contents to get a stable result
        std::mt19937 rng(dev());
        std::uniform_real_distribution<double> dist(0,1);
        for (uint64_t i = 0; i < 2*n; ++i) {
            X[i] = dist(rng);
        }
        // do layout
        sgd2::layout_sparse_unweighted(n, X.data(), I.size(), I.data(), J.data(), pivots, t_max, eps);
        for (uint64_t i = 0; i < 2*n; i+=2) {
            //std::cerr << "i = " << i << std::endl;
            uint64_t j = component_ids[i/2]-1;
            layout[j] = X[i] + max_x;
            layout[j+1] = X[i+1];
        }
        // set new max_x
        for (uint64_t i = 0; i < 2*n; i+=2) {
            max_x = std::max(X[i], max_x);
        }
        max_x += x_padding;
    }
    // all the weakly connected component layouts have been merged here
    return layout;
}

}
}
