#ifndef dgraph_dynamic_types_hpp
#define dgraph_dynamic_types_hpp

#include "dynamic.hpp"

namespace odgi {

typedef dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,256,16> > suc_bv;
typedef dyn::lciv<dyn::packed_vector,512,1> lciv_iv;
typedef dyn::wt_string<dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,256,16> > > wt_str;
typedef dyn::hacked_vector suc_iv;

}

#endif
