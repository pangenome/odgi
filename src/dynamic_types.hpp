#ifndef dgraph_dynamic_types_hpp
#define dgraph_dynamic_types_hpp

#include "dynamic.hpp"

namespace dsgvg {

typedef dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,256,16> > suc_bv;
//typedef dyn::spsi<dyn::packed_vector,1024,1> spsi_iv;
typedef dyn::lciv<dyn::packed_vector,256,16> lciv_iv;
typedef dyn::wt_string<dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,256,16> > > wt_str;
typedef dyn::hacked_vector suc_iv;

}

#endif
