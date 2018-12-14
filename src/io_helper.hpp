#ifndef gfa_io_helper
#define gfa_io_helper
#include "gfakluge.hpp"
#include "graph.hpp"
#include "handle_types.hpp"


inline void dank_to_gfa_stream(const graph_t& sd, std::ostream& os) const {
    
};

inline graph_t dank_from_gfa_file(const std::string& filename) const{

    gfak::GFAKluge gg;
    gg.parse_gfa_file(std::string(filename));

};



#endif
