#include "crush_n.hpp"

namespace odgi {
namespace algorithms {

void crush_n(odgi::graph_t& graph) {
    graph.for_each_handle([&](const handle_t& handle) {
        // strip Ns from start
        const std::string src = graph.get_sequence(handle);
        std::string seq;
        seq.reserve(src.size());
        bool in_n = false;
        for (auto c : src) {
            if (c == 'N') {
                if (in_n) {
                    continue;
                } else {
                    in_n = true;
                }
            } else {
                in_n = false;
            }
            seq.push_back(c);
        }
        graph.set_handle_sequence(handle, seq);
    }, true); // in parallel
}

}
}
