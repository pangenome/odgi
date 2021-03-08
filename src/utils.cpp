#include <string>
#include <algorithm>
#include "utils.hpp"

namespace utils {
    bool is_number(const std::string &s) {
        return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
    }
}