#ifndef DSGVG_HASH_HPP_INCLUDED
#define DSGVG_HASH_HPP_INCLUDED

#include <cstdint>

namespace odgi {

uint32_t djb2_hash32(const char *str)
{
    uint32_t hash = 5381;
    int c;

    while (c = *str++)
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}

uint64_t djb2_hash64(const char *str)
{
    uint64_t hash = 5381;
    int c;

    while (c = *str++)
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}

}

#endif
