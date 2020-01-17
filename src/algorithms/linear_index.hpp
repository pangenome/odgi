#pragma once

#include <vector>
#include <string>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <cassert>

namespace odgi {
namespace algorithms {

using namespace handlegraph;

class linear_index_t {
public:
    std::string graph_seq;
    std::vector<uint64_t> handle_positions;
    uint64_t position_of_handle(const handle_t& handle);
    linear_index_t(const PathHandleGraph& graph);
};

}
}
