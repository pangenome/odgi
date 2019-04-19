#include "coverage.hpp"

namespace odgi {
namespace algorithms {

std::vector<handle_t> find_handles_exceeding_coverage_limits(const MutablePathDeletableHandleGraph& graph, uint64_t min_coverage, uint64_t max_coverage) {
    std::vector<handle_t> handles;
    graph.for_each_handle([&](const handle_t& handle) {
            uint64_t step_count = graph.occurrences_of_handle(handle).size();
            if (min_coverage && step_count < min_coverage || max_coverage && step_count > max_coverage) {
                handles.push_back(handle);
            }
        });
    return handles;
}

}
}
