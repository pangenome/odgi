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
    // translate the order
    std::reverse(order.begin(), order.end());
    ska::flat_hash_set<handlegraph::nid_t> seen;
    std::vector<handle_t> translated_order;
    for (auto& handle : order) {
        //assert(dagified_to_orig.find(into.get_id(handle)) != dagified_to_orig.end());
        auto vs_orig = dagified_to_orig(into.get_id(handle));
        handlegraph::nid_t id = vs_orig.first;
        //assert(id > 0);
        if (!seen.count(id)) {
            translated_order.push_back(base.get_handle(id));
            seen.insert(id);
        }
    }
    std::reverse(translated_order.begin(), translated_order.end());
    assert(translated_order.size() == base.get_node_count());
    return translated_order;
}

}
}
