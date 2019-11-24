#include "dagify_sort.hpp"

namespace odgi {
namespace algorithms {

std::vector<handle_t> dagify_sort(const HandleGraph& base, MutableHandleGraph& split, MutableHandleGraph& into) {
    auto split_to_orig = algorithms::split_strands(&base, &split);
    auto dagified_to_split = algorithms::dagify(&split, &into, 1);
    auto dagified_to_orig = [&](handlegraph::nid_t id) {
        return split_to_orig[dagified_to_split[id]];
    };
    auto order = algorithms::topological_order(&into, true, false);
    // find the mean position in the order for each original handle
    ska::flat_hash_map<handlegraph::nid_t, std::pair<uint64_t, uint64_t>> pos_accum;
    for (uint64_t i = 0; i < order.size(); ++i) {
        auto& handle = order[i];
        auto e = dagified_to_orig(into.get_id(handle));
        if (e.second) continue;
        handlegraph::nid_t id = e.first;
        pos_accum[id].first += i;
        ++pos_accum[id].second;
    }
    // sort the original ids by their average position in the dagified sort
    std::vector<std::pair<handlegraph::nid_t, double>> avg_pos; avg_pos.reserve(pos_accum.size());
    for (auto& e : pos_accum) {
        avg_pos.push_back(std::make_pair(e.first, (double)e.second.first/(double)e.second.second));
    }
    std::sort(avg_pos.begin(), avg_pos.end(),
              [](const std::pair<handlegraph::nid_t, double>& a,
                 const std::pair<handlegraph::nid_t, double>& b) {
                  return a.second < b.second;
              });
    std::vector<handle_t> translated_order;
    for (auto& e : avg_pos) {
        translated_order.push_back(base.get_handle(e.first));
    }
    //std::reverse(translated_order.begin(), translated_order.end());
    assert(translated_order.size() == base.get_node_count());
    return translated_order;
}

}
}
