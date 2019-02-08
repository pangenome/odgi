#ifndef dgraph_handle_helper
#define dgraph_handle_helper
#include <cstdint>
#include <cassert>
#include "handle_types.hpp"

namespace dsgvg {
/// View a handle as an integer
inline uint64_t& as_integer(handle_t& handle) {
    return reinterpret_cast<uint64_t&>(handle);
}

/// View a const handle as a const integer
inline const uint64_t& as_integer(const handle_t& handle) {
    return reinterpret_cast<const uint64_t&>(handle);
}

/// View an integer as a handle
inline handle_t& as_handle(uint64_t& value) {
    return reinterpret_cast<handle_t&>(value);
}

/// View a const integer as a const handle
inline const handle_t& as_handle(const uint64_t& value) {
    return reinterpret_cast<const handle_t&>(value);
}

/// Define equality on handles
inline bool operator==(const handle_t& a, const handle_t& b) {
    return as_integer(a) == as_integer(b);
}

/// Define inequality on handles
inline bool operator!=(const handle_t& a, const handle_t& b) {
    return as_integer(a) != as_integer(b);
}

struct handle_helper{

    /// Extract the packed integer
    inline static uint64_t unpack_number(const handle_t& handle) {
        return as_integer(handle) >> 1;
    }
    
    /// Extract the packed bit
    inline static bool unpack_bit(const handle_t& handle) {
        return as_integer(handle) & 1;
    }
    
    /// Pack up an integer and a bit into a handle
    inline static handle_t pack(const uint64_t& number, const bool& bit) {
        // Make sure the number doesn't use all the bits
        assert(number < (0x1ULL << 63));
        
        return as_handle((number << 1) | (bit ? 1 : 0));
    }
    
    /// Toggle the packed bit and return a new handle
    inline static handle_t toggle_bit(const handle_t& handle) {
        return as_handle(as_integer(handle) ^ 1);
    }
};

///
/// Path handles
///

/// View a path handle as an integer
inline uint64_t& as_integer(path_handle_t& handle) {
    return reinterpret_cast<uint64_t&>(handle);
}

/// View a const path handle as a const integer
inline const uint64_t& as_integer(const path_handle_t& handle) {
    return reinterpret_cast<const uint64_t&>(handle);
}

/// View an integer as a path handle
inline path_handle_t& as_path_handle(uint64_t& value) {
    return reinterpret_cast<path_handle_t&>(value);
}

/// View a const integer as a const path handle
inline const path_handle_t& as_path_handle(const uint64_t& value) {
    return reinterpret_cast<const path_handle_t&>(value);
}

/// Define equality on path handles
inline bool operator==(const path_handle_t& a, const path_handle_t& b) {
    return as_integer(a) == as_integer(b);
}

/// Define inequality on path handles
inline bool operator!=(const path_handle_t& a, const path_handle_t& b) {
    return as_integer(a) != as_integer(b);
}

///
/// Occurrence handles
///

/// View an occurrence handle as an integer
inline int64_t* as_integers(occurrence_handle_t& occurrence_handle) {
    return reinterpret_cast<int64_t*>(&occurrence_handle);
}

/// View a const occurrence handle as a const integer
inline const int64_t* as_integers(const occurrence_handle_t& occurrence_handle) {
    return reinterpret_cast<const int64_t*>(&occurrence_handle);
}

/// Define equality on occurrence handles
inline bool operator==(const occurrence_handle_t& a, const occurrence_handle_t& b) {
    return as_integers(a)[0] == as_integers(b)[0] && as_integers(a)[1] == as_integers(b)[1];
}

/// Define inequality on occurrence handles
inline bool operator!=(const occurrence_handle_t& a, const occurrence_handle_t& b) {
    return !(a == b);
}

}


#endif
