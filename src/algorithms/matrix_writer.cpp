#include "matrix_writer.hpp"

namespace odgi {
namespace algorithms {

void write_as_sparse_matrix(std::ostream& out, const PathHandleGraph& graph, bool weight_by_edge_depth, bool weight_by_edge_delta) {
    uint64_t edge_count = 0;
    graph.for_each_edge([&](const edge_t& edge) {
            ++edge_count;
        });
    out << graph.max_node_id() << " " << graph.max_node_id() << " " << edge_count*2 << std::endl;
    graph.for_each_edge([&](const edge_t& edge) {
            // how many paths cross the edge?
            double weight = 0;
            if (weight_by_edge_depth) {
                graph.for_each_step_on_handle(edge.first, [&](const step_handle_t& step) {
                        if (graph.get_handle_of_step(step) == edge.first && graph.has_next_step(step)) {
                            const handle_t next = graph.get_handle_of_step(graph.get_next_step(step));
                            if (next == edge.second) {
                                ++weight;
                            }
                        } else if (graph.get_handle_of_step(step) == graph.flip(edge.first) && graph.has_previous_step(step)) {
                            const handle_t prev = graph.get_handle_of_step(graph.get_previous_step(step));
                            if (prev == graph.flip(edge.second)) {
                                ++weight;
                            }
                        }
                    });
            } else {
                weight = 1;
            }
            if (weight_by_edge_delta) {
                double delta = std::abs(graph.get_id(edge.first) - graph.get_id(edge.second));
                if (delta == 0) delta = 1;
                weight = 1 / delta;
            }
            out << graph.get_id(edge.first) << " " << graph.get_id(edge.second) << " " << weight << std::endl;
            out << graph.get_id(edge.second) << " " << graph.get_id(edge.first) << " " << weight << std::endl;
        });
}

}
}
