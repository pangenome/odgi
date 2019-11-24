#include "random_order.hpp"

namespace odgi {

namespace algorithms {

std::vector<handle_t> random_order(const HandleGraph& graph) {
    std::vector<handle_t> order;
    graph.for_each_handle([&order](const handle_t& handle) {
            order.push_back(handle);
        });
    std::random_device dev;
    std::mt19937 rng(dev());
    std::shuffle(order.begin(), order.end(), rng);
    return order;
}

}
}
