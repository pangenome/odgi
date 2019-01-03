#ifndef dgraph_dynamic_types_hpp
#define dgraph_dynamic_types_hpp

#include "dynamic.hpp"

namespace dankgraph {

typedef dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,8192,16> > suc_bv;
//typedef dyn::gap_bitvector<dyn::spsi<dyn::packed_vector,8192,16> > suc_bv;
typedef dyn::wt_string<suc_bv> wt_str;
//typedef dyn::wt_string<dyn::rle_str> wt_str;

}

#endif
