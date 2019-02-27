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

// Utility functions for making the compiler do what we want.

#include <cstdint>
#include <cstring>
#include <vector>

namespace varint {

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

inline std::vector<uint8_t> encode(const std::vector<uint64_t> &in) {
    std::vector<uint8_t> out;
    return encode(in, out);
}

inline void decode(const uint8_t *in, uint64_t *out, size_t count) {
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
}

inline uint64_t length(const uint64_t& i) {
    if (i < cut1) {
        return 1;
    } else if (i <= cut1 + 255 + 256 * (cut2 - 1 - cut1)) {
        return 2;
    } else {
        return (64 - count_leading_zeros_64(i) + 7) / 8;
    }
}

}
