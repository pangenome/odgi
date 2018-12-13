#ifndef dank_path_hpp
#define dank_path_hpp

#include "dynamic.hpp"

namespace dankgraph {

struct step_t {
    uint64_t start = 0;
    uint64_t length = 0;
    bool strand = 0;
};

class path_t {
public:
    std::string name;
    path_t(const std::string& n) { name = n; }
    void append_occurrence(const uint64_t& start, const uint64_t& length, bool strand);
    void clear(void);
    uint64_t occurrence_count(void) const;
    step_t get_occurrence(uint64_t step) const;
private:
    dyn::wt_string<dyn::suc_bv> starts;
    dyn::wt_string<dyn::suc_bv> lengths;
    dyn::suc_bv strands;
};

inline void path_t::clear(void) {
    dyn::wt_string<dyn::suc_bv> null_wt;
    dyn::suc_bv null_bv;
    starts = null_wt;
    lengths = null_wt;
    strands = null_bv;
}

inline void path_t::append_occurrence(const uint64_t& start,
                                      const uint64_t& length,
                                      bool strand) {
    starts.push_back(start);
    lengths.push_back(length);
    strands.push_back(strand);
}

inline uint64_t path_t::occurrence_count(void) const {
    return starts.size();
}

inline step_t path_t::get_occurrence(uint64_t step) const {
    step_t result;
    result.start = starts.at(step);
    result.length = lengths.at(step);
    result.strand = strands.at(step);
    return result;
}

}

#endif
