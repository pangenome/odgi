#ifndef dank_path_hpp
#define dank_path_hpp

#include "dynamic.hpp"

namespace dankgraph {

class path_t {
public:
    std::string name;
    path_t(const std::string& n) { name = n; }
    void append_occurrence(const uint64_t& start, const uint64_t& length, bool strand);
    void clear(void);
private:
    dyn::wt_string<dyn::suc_bv> starts;
    dyn::wt_string<dyn::suc_bv> lengths;
    dyn::suc_bv strands;
};

inline void path_t::append_occurrence(const uint64_t& start,
                                      const uint64_t& length,
                                      bool strand) {
    starts.push_back(start);
    lengths.push_back(length);
    strands.push_back(strand);
}

}

#endif
