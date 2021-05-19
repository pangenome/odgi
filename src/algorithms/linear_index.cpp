#include "linear_index.hpp"

namespace odgi {
namespace algorithms {

linear_index_t::linear_index_t(const PathHandleGraph& graph) {
    // generate our flattened sequence vector and a positional mapping into it for each handle
    uint64_t graph_seq_size = 0;
    graph.for_each_handle([&](const handle_t& h) {
                              graph_seq_size += graph.get_length(h);
                          });
    graph_seq.reserve(graph_seq_size);
    handle_positions.reserve(graph.get_node_count());
    uint64_t curr_pos_in_seq = 0;
    graph.for_each_handle([&](const handle_t& h) {
        graph_seq.append(graph.get_sequence(h));
        // verify that our graph handle space is compact
        // it should be when using a freshly loaded odgi graph
        assert(number_bool_packing::unpack_number(h) == handle_positions.size());
        handle_positions.push_back(curr_pos_in_seq);
        curr_pos_in_seq += graph.get_length(h);
    });
}

uint64_t linear_index_t::position_of_handle(const handle_t& handle) {
    return handle_positions.at(number_bool_packing::unpack_number(handle));
}

}
}
