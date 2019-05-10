#include "id_ordered_paths.hpp"

namespace odgi {
namespace algorithms {

std::vector<path_handle_t> id_ordered_paths(const PathHandleGraph& g, bool rev) {
    // accumulate the paths by their minimum id
    std::vector<std::pair<uint64_t, path_handle_t> > paths;
    g.for_each_path_handle([&](const path_handle_t& p) {
            // get the min id
            uint64_t min_id = std::numeric_limits<uint64_t>::max();
            g.for_each_step_in_path(p, [&](const step_handle_t& occ) {
                    min_id = std::min((uint64_t)g.get_id(g.get_handle_of_step(occ)), min_id);
                });
            paths.push_back(make_pair(min_id, p));
        });
    std::sort(paths.begin(), paths.end(), [&](const std::pair<uint64_t, path_handle_t>& a,
                                              const std::pair<uint64_t, path_handle_t>& b) {
                  return a.first < b.first || a.first == b.first && as_integer(a.second) < as_integer(b.second);
              });
    if (rev) {
        std::reverse(paths.begin(), paths.end());
    }
    std::vector<path_handle_t> order; order.reserve(paths.size());
    for (auto& p : paths) {
        order.push_back(p.second);
    }
    return order;
}

}
}
