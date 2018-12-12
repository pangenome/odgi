#ifndef gfa_io_helper
#define gfa_io_helper
#include "gfakluge.hpp"
#include "handle_types.hpp"


inline void dank_to_gfa_stream(SuccinctDynamicSequenceGraph* sd, std::ostream& os) const {
    
};

inline void dank_from_gfa_file(char* filename, SuccinctDynamicSequenceGraph* sd) const{

    gfak::GFAKluge gg;
    gg.parse_gfa_file(std::string(filename));

};



#endif