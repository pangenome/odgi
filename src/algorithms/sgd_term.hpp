#pragma once

#include <handlegraph/handle_graph.hpp>

namespace odgi {
    namespace algorithms {

        using namespace handlegraph;

        struct sgd_term_t {
            handle_t i, j;
            double d, w;

            sgd_term_t(const handle_t &i, const handle_t &j, const double &d, const double &w) : i(i), j(j), d(d),
                                                                                                 w(w) {}
        };
    }
}