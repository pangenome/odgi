#include "odgi.hpp"
#include "gfa_to_handle.hpp"

#include <filesystem>

using namespace handlegraph;

#ifndef ODGI_UTILS_H
#define ODGI_UTILS_H

namespace utils {
    bool is_number(const std::string &s);
    std::string to_string_custom(double value);
    void graph_deep_copy(const odgi::graph_t &source,
                         odgi::graph_t* target);
	bool ends_with(const std::string &fullString, const std::string &ending);
	int handle_gfa_odgi_input(const std::string infile, const std::string subcommmand_name, const bool progress,
							  const uint64_t num_threads, odgi::graph_t &graph);
	/// this function will return n % d
	/// it is assumed that d is one of 1, 2, 4, 8, 16, 32, ....
	uint64_t modulo(const uint64_t n, const uint64_t d);
}


#endif //ODGI_UTILS_H
