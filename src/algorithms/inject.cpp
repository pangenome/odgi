#include "inject.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

void inject_ranges(MutablePathHandleGraph& graph,
                   const ska::flat_hash_map<path_handle_t, std::pair<std::string, std::vector<interval_t>>>& path_intervals) {
    // first we collect cut points
    ska::flat_hash_map<handle_t, std::vector<size_t>> cut_points;
    // then we cut the nodes in the graph
    // then we iterate back through the sorted path intervals and add paths at the appropriate points
}

void chop_at(MutablePathDeletableHandleGraph &graph,
             const ska::flat_hash_map<handle_t, std::vector<size_t>>& cut_points,
             const uint64_t &nthreads, const bool &show_info) {

    std::vector<std::tuple<uint64_t, uint64_t, handle_t>> originalRank_inChoppedNodeRank_handle;
    std::vector<std::pair<uint64_t, handle_t>> originalRank_handleToChop;

    uint64_t rank = 0;
    graph.for_each_handle([&](const handle_t &handle) {
        if (cut_points.find(handle) != cut_points.end()) {
            originalRank_handleToChop.push_back(std::make_pair(rank, handle));
        } else {
            originalRank_inChoppedNodeRank_handle.push_back(std::make_tuple(rank, 0, handle));
        }
        rank++;
    });

    if (show_info) {
        std::cerr << "[odgi::inject::chop_at] " << originalRank_handleToChop.size() << " node(s) to chop." << std::endl;
    }

    for (auto rank_handle : originalRank_handleToChop) {
        // get divide points
        auto& handle = rank_handle.second;
        uint64_t length = graph.get_length(handle);
        auto f = cut_points.find(handle);
        assert(f != cut_points.end());
        const std::vector<size_t>& cut_offsets = f->second;
        rank = 0;
        for (auto chopped_handle : graph.divide_handle(rank_handle.second, cut_offsets)) {
            originalRank_inChoppedNodeRank_handle.push_back(
                std::make_tuple(rank_handle.first, rank, chopped_handle));
            rank++;
        }
    }

    ips4o::parallel::sort(
        originalRank_inChoppedNodeRank_handle.begin(), originalRank_inChoppedNodeRank_handle.end(),
        std::less<>(), nthreads
        );

    std::vector<handle_t> new_handles;
    for (auto x_y_z : originalRank_inChoppedNodeRank_handle) {
        new_handles.push_back(std::get<2>(x_y_z));
    }

    graph.apply_ordering(new_handles, true);
}


}

}
