#include <odgi.hpp>

using namespace handlegraph;

#ifndef ODGI_UTILS_H
#define ODGI_UTILS_H

namespace utils {
    bool is_number(const std::string &s);
    void graph_deep_copy(const odgi::graph_t &source,
                         odgi::graph_t* target);
}


#endif //ODGI_UTILS_H
