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
                // might be more correct if this were xor
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

uint64_t cut_tips(
    MutablePathDeletableHandleGraph& graph,
    uint64_t min_coverage) {
    auto tips = tip_handles(graph);
    std::vector<handle_t> drop_tips;
    std::sort(tips.begin(), tips.end());
    for (auto& tip : tips) {
        std::vector<step_handle_t> to_destroy;
        graph.for_each_step_on_handle(
            tip,
            [&](const step_handle_t& step) {
                to_destroy.push_back(step);
            });
        if (!min_coverage || to_destroy.size() < min_coverage) {
            // odgi quirk / baddness: we have to destroy step handles from largest to smallest
            // a generic implementation would iterate through step handles while erasing
            std::sort(to_destroy.begin(), to_destroy.end(),
                      [](const step_handle_t& a, const step_handle_t& b) {
                          return as_integers(a)[1] > as_integers(b)[1]; });
            for (auto& step : to_destroy) {
                graph.rewrite_segment(step, step, {});
            }
            drop_tips.push_back(tip);
        }
    }
    for (auto& tip : drop_tips) {
        graph.destroy_handle(tip);
    }
    return tips.size();
}

}

}
