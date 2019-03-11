#include <cassert>
#include "varint.hpp"
#include <iostream>

namespace odgi {

namespace varint {

static const uint8_t MSB = 0x80;
static const uint8_t MSBALL = ~0x7F;

static const uint64_t N1 = 128; // 2 ^ 7
static const uint64_t N2 = 16384;
static const uint64_t N3 = 2097152;
static const uint64_t N4 = 268435456;
static const uint64_t N5 = 34359738368;
static const uint64_t N6 = 4398046511104;
static const uint64_t N7 = 562949953421312;
static const uint64_t N8 = 72057594037927936;
static const uint64_t N9 = 9223372036854775808U;

uint8_t length(uint64_t n) {
    return (
        n < N1 ? 1
        : n < N2 ? 2
        : n < N3 ? 3
        : n < N4 ? 4
        : n < N5 ? 5
        : n < N6 ? 6
        : n < N7 ? 7
        : n < N8 ? 8
        : n < N9 ? 9
        :         10
        );
}

uint64_t length(const uint64_t* v, uint64_t n) {
    uint64_t len = 0;
    for ( ; n > 0; --n) {
        len += length(*v++);
    }
    return len;
}

uint64_t length(const std::vector<uint64_t>& v) {
    return length((uint64_t*)v.data(), v.size());
}

uint8_t* encode(uint64_t n, uint8_t* ptr) {
    while (n & MSBALL) {
        *(ptr++) = (n & 0xFF) | MSB;
        n = n >> 7;
    }
    *ptr++ = n;
    return ptr;
}

uint8_t* encode(const uint64_t* in, uint8_t* ptr, uint64_t c) {
    for ( ; c > 0; --c) {
        ptr = encode(*in++, ptr);
    }
    return ptr;
}

uint8_t* encode(const std::vector<uint64_t>& v, uint8_t* ptr) {
    return encode((uint64_t*)v.data(), ptr, v.size());
}

uint8_t* decode(uint64_t* out, uint8_t* ptr) {
    uint8_t bits = 0;
    uint64_t ll = 0;
    while (*ptr & MSB) {
        ll = *ptr;
        *out += ((ll & 0x7F) << bits);
        ptr++;
        bits += 7;
    }
    ll = *ptr++;
    *out += ((ll & 0x7F) << bits);
    return ptr;
}

uint8_t* decode(uint64_t* out, uint8_t* ptr, uint64_t c) {
    for ( ; c > 0; --c) {
        ptr = decode(out++, ptr);
    }
    return ptr;
}

uint8_t* seek(uint8_t* ptr, uint64_t n) {
    for ( ; n > 0; --n) {
        while (*ptr & MSB) {
            ptr++;
        }
        ptr++;
    }
    return ptr;
}

uint64_t bytes(uint8_t* ptr, uint64_t n) {
    return seek(ptr, n) - ptr;
}

}

}
