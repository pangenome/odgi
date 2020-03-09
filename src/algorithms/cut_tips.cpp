#include "cut_tips.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

std::vector<handle_t> tip_handles(
    const HandleGraph& graph) {
    std::vector<handle_t> tips;
    graph.for_each_handle(
        [&tips,&graph](const handle_t& handle) {
            // check if one end of the handle is a tip
            if (graph.get_degree(handle, false) == 0
                || graph.get_degree(handle, true) == 0) {
                tips.push_back(handle);
            }
        });
    return tips;
}

uint64_t cut_tips(
    DeletableHandleGraph& graph) {
    auto tips = tip_handles(graph);
    for (auto& tip : tips) {
        graph.destroy_handle(tip);
    }
    return tips.size();
}

}

}
