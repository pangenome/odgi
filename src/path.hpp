#ifndef dank_path_hpp
#define dank_path_hpp

#include "dynamic_structs.hpp"
#include "sdsl/bit_vectors.hpp"

namespace dankgraph {

class path_t {
public:
    std::string name;
    path_t(const std::string& n) { name = n; }
    void append_occurrence(const uint64_t& start, const uint64_t& length, bool strand);
    void clear(void);
private:
    //size_t filled = 0;
    SuccinctDynamicVector starts;
    SuccinctDynamicVector lengths;
    SuccinctDynamicVector strands;
};

inline void path_t::append_occurrence(const uint64_t& start,
                                      const uint64_t& length,
                                      bool strand) {
    starts.append(start);
    lengths.append(length);
    strands.append(strand);
}

inline void path_t::clear(void) {
    starts.clear();
    lengths.clear();
    strands.clear();
}

}

#endif
