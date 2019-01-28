#ifndef dgraph_dynamic_types_hpp
#define dgraph_dynamic_types_hpp

#include "dynamic.hpp"

namespace dankgraph {

typedef dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,2048,2> > suc_bv;
//typedef dyn::spsi<dyn::packed_vector,1024,1> spsi_iv;
typedef dyn::lciv<dyn::packed_vector,1024,1> lciv_iv;
typedef dyn::wt_string<dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,2048,2> > > wt_str;

}

#endif
