#include "diffpriv.hpp"

namespace odgi {
namespace algorithms {

void diff_priv(
    const PathHandleGraph& graph,
    PathHandleGraph& priv,
    double epsilon,
    double target_coverage,
    uint64_t bp_limit) {

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<uint64_t> dist(1, graph.get_node_count());
    std::uniform_int_distribution<uint64_t> coin(0, 1);

    typedef std::vector<std::pair<step_handle_t, step_handle_t>> step_ranges_t;

    // algorithm
    while (true) {
        // we randomly sample a starting node (todo: step)
        uint64_t rand_id = dist(mt);
        // we collect all steps on the node, picking a random orientation
        handle_t h = graph.get_handle(rand_id, (bool)coin(mt));
        // we collect all potential forward extensions
        step_ranges_t ranges;
        graph.for_each_step_on_handle(
            h, [&](const step_handle_t& s) {
                ranges.push_back(std::make_pair(s, s));
            });
        // sampling loop
        while (!ranges.empty()) {
            // next handles
            std::map<handle_t, step_ranges_t> nexts;
            for (auto& range : ranges) {
                auto& s = range.second;
                if (graph.has_next_step(s)) {
                    step_handle_t q = graph.get_next_step(s);
                    handle_t n = graph.get_handle_of_step(q);
                    nexts[n].push_back(std::make_pair(range.first, q));
                }
            }
            // we calculate the utilities of each potential next
            // and then apply the exponential mechanism to sample
        }
    }

    // and look for potential extensions of open path intervals

}

}
}
