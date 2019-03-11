#ifndef _ODGI_VARINT_H_
#define _ODGI_VARINT_H_

#include <cstdint>
#include <cstring>
#include <vector>

namespace odgi {

// Copyright 2016 Jakob Stoklund Olesen
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

namespace sqvarint {

// Utility functions for making the compiler do what we want.

// Count the number of contiguous zero bits starting from the MSB.
// It is undefined behavior if x is 0.
inline unsigned
count_leading_zeros_64(uint64_t x)
{
  return __builtin_clzll(x);
}

// Count the number of contiguous zero bits starting from the LSB.
// It is undefined behavior if x is 0.
inline unsigned
count_trailing_zeros_32(uint32_t x)
{
  return __builtin_ctz(x);
}

// Load an unsigned 64-bit number from an unaligned address.
// This will usually be translated into one instruction.
inline uint64_t
unaligned_load_u64(const uint8_t* p)
{
  uint64_t x;
  std::memcpy(&x, p, 8);
  return x;
}

inline uint64_t
unaligned_load(const uint8_t* p, uint8_t s)
{
    uint64_t x;
    std::memcpy(&x, p, s);
    return x;
}

// Load an unsigned 16-bit number from an unaligned address.
// This will usually be translated into one instruction.
inline uint16_t
unaligned_load_u16(const uint8_t* p)
{
  uint16_t x;
  std::memcpy(&x, p, 2);
  return x;
}

// Hint to the compiler that X is probably true.
#define LIKELY(X) __builtin_expect(!!(X), 1)

// This implements a variant of the SQLite variable-length integer encoding.
//
// The first byte B0 determines the length:
//
// 0-184: 1 byte, value = B0.
// 185-248: 2 bytes, value = 185 + B1 + 256*(B0-185).
// 249-255: 3-9 bytes, B0-249+2 little-endian encoded bytes following B0.

//#include "varint.h"

const unsigned cut1 = 185;
const unsigned cut2 = 249;

inline std::vector<uint8_t> encode(const std::vector<uint64_t> &in, std::vector<uint8_t>& out) {
    for (auto x : in) {
        if (x < cut1) {
            // 1 byte.
            out.push_back(x);
        } else if (x <= cut1 + 255 + 256 * (cut2 - 1 - cut1)) {
            // 2 bytes.
            x -= cut1;
            out.push_back(cut1 + (x >> 8));
            out.push_back(x & 0xff);
        } else {
            // 3-9 bytes.
            unsigned bits = 64 - count_leading_zeros_64(x);
            unsigned bytes = (bits + 7) / 8;
            out.push_back(cut2 + (bytes - 2));
            for (unsigned n = 0; n < bytes; n++) {
                out.push_back(x & 0xff);
                x >>= 8;
            }
        }
    }
    return out;
}

inline uint8_t* encode(const std::vector<uint64_t> &in, uint8_t* out) {
    for (auto x : in) {
        if (x < cut1) {
            // 1 byte.
            *out++ = x;
        } else if (x <= cut1 + 255 + 256 * (cut2 - 1 - cut1)) {
            // 2 bytes.
            x -= cut1;
            *out++ = cut1 + (x >> 8);
            *out++ = x & 0xff;
        } else {
            // 3-9 bytes.
            unsigned bits = 64 - count_leading_zeros_64(x);
            unsigned bytes = (bits + 7) / 8;
            *out++ = cut2 + (bytes - 2);
            for (unsigned n = 0; n < bytes; n++) {
                *out++ = x & 0xff;
                x >>= 8;
            }
        }
    }
    return out;
}

inline std::vector<uint8_t> encode(const std::vector<uint64_t> &in) {
    std::vector<uint8_t> out;
    return encode(in, out);
}

inline uint8_t* encode(const uint64_t* in, uint8_t* out, uint64_t count) {
    for (uint64_t i = 0; i < count; ++i) {
        uint64_t x = *in++;
        if (x < cut1) {
            // 1 byte.
            *out++ = x;
        } else if (x <= cut1 + 255 + 256 * (cut2 - 1 - cut1)) {
            // 2 bytes.
            x -= cut1;
            *out++ = cut1 + (x >> 8);
            *out++ = x & 0xff;
        } else {
            // 3-9 bytes.
            unsigned bits = 64 - count_leading_zeros_64(x);
            unsigned bytes = (bits + 7) / 8;
            *out++ = cut2 + (bytes - 2);
            for (unsigned n = 0; n < bytes; n++) {
                *out++ = x & 0xff;
                x >>= 8;
            }
        }
    }
    return out;
}

inline uint8_t* decode(uint64_t *out, uint8_t *in, size_t count) {
    while (count-- > 0) {
        uint8_t b0 = *in++;
        if (LIKELY(b0 < cut1)) {
            *out++ = b0;
        } else if (b0 < cut2) {
            uint8_t b1 = *in++;
            *out++ = cut1 + b1 + ((b0 - cut1) << 8);
        } else {
            size_t sh = b0 - cut2;
            *out++ = unaligned_load_u64(in) & ((uint64_t(1) << 8 * sh << 16) - 1);
            in += 2 + sh;
        }
    }
    return in;
}

inline uint64_t length(const uint64_t& i) {
    if (i < cut1) {
        return 1;
    } else if (i <= cut1 + 255 + 256 * (cut2 - 1 - cut1)) {
        return 2;
    } else {
        return 1 + (64 - count_leading_zeros_64(i) + 7) / 8;
    }
}

inline uint64_t length(const std::vector<uint64_t>& v) {
    uint64_t l = 0;
    for (auto& i : v) l += length(i);
    return l;
}

inline uint64_t length(const uint64_t* v, uint64_t count) {
    uint64_t l = 0;
    for (uint64_t i = 0; i < count; ++i) l += length(*v++);
    return l;
}

inline uint8_t* seek(uint8_t *in, size_t count) {
    while (count-- > 0) {
        uint8_t b0 = *in++;
        if (LIKELY(b0 < cut1)) {
        } else if (b0 < cut2) {
            uint8_t b1 = *in++;
        } else {
            size_t sh = b0 - cut2;
            in += 2 + sh;
        }
    }
    return in;
}

inline uint64_t bytes(uint8_t *in, size_t count) {
    uint64_t bytes = 0;
    while (count-- > 0) {
        uint8_t b0 = *in++;
        bytes++;
        if (LIKELY(b0 < cut1)) {
        } else if (b0 < cut2) {
            uint8_t b1 = *in++;
            bytes++;
        } else {
            size_t sh = b0 - cut2;
            in += 2 + sh;
            bytes += 2 + sh;
        }
    }
    return bytes;
}

}

namespace msbvarint {

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

inline uint64_t length(uint64_t n) {
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

inline uint64_t length(const uint64_t* v, uint64_t n) {
    uint64_t len = 0;
    for ( ; n > 0; --n) {
        len += length(*v++);
    }
    return len;
}

inline uint64_t length(const std::vector<uint64_t>& v) {
    return length((uint64_t*)v.data(), v.size());
}

inline uint8_t* encode(uint64_t n, uint8_t* ptr) {
    uint8_t* buf = ptr;
    uint64_t l = length(n);
    for ( ; l > 1; --l) {
        // why this broken?
        //while (n & MSBALL) {
        *(ptr++) = (n & 0xFF) | MSB;
        n = n >> 7;
    }
    *ptr++ = n;
    return ptr;
}

inline uint8_t* encode(const uint64_t* in, uint8_t* ptr, uint64_t c) {
    for ( ; c > 0; --c) {
        ptr = encode(*in++, ptr);
    }
    return ptr;
}

inline uint8_t* encode(const std::vector<uint64_t>& v, uint8_t* ptr) {
    return encode((uint64_t*)v.data(), ptr, v.size());
}

inline uint8_t* decode(uint64_t* out, uint8_t* ptr) {
    uint8_t bits = 0;
    uint64_t ll = 0;
    while (*ptr & MSB) {
        ll = *ptr++;
        *out += ((ll & 0x7F) << bits);
        bits += 7;
    }
    ll = *ptr++;
    *out += ((ll & 0x7F) << bits);
    return ptr;
}

inline uint8_t* decode(uint64_t* out, uint8_t* ptr, uint64_t c) {
    for ( ; c > 0; --c) {
        ptr = decode(out++, ptr);
    }
    return ptr;
}

inline uint8_t* seek(uint8_t* ptr, uint64_t n) {
    for ( ; n > 0; --n) {
        while (*ptr & MSB) {
            ptr++;
        }
        ptr++;
    }
    return ptr;
}

inline uint64_t bytes(uint8_t* ptr, uint64_t n) {
    return seek(ptr, n) - ptr;
}

}

}

#endif
