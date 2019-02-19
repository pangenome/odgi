#ifndef DSGVG_PHF_HPP
#define DSGVG_PHF_HPP

#include <cstdint>
#include "BooPHF.h"

namespace dsgvg {

typedef boomphf::SingleHashFunctor<uint64_t>  hasher_t;
typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;

}

#endif
