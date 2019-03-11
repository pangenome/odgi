#ifndef _ODGI_VARINT_H_
#define _ODGI_VARINT_H_

#include <cstdint>
#include <vector>

namespace odgi {

namespace varint {

uint64_t length(uint64_t n);
uint64_t length(const uint64_t* v, uint64_t n);
uint64_t length(const std::vector<uint64_t>& v);
uint8_t* encode(uint64_t n, uint8_t* ptr);
uint8_t* encode(const uint64_t* in, uint8_t* ptr, uint64_t c);
uint8_t* encode(const std::vector<uint64_t>& v, uint8_t* ptr);
uint8_t* decode(uint64_t* out, uint8_t* ptr);
uint8_t* decode(uint64_t* out, uint8_t* ptr, uint64_t c);
uint8_t* seek(uint8_t* ptr, uint64_t n);
uint64_t bytes(uint8_t* ptr, uint64_t n);

}

}

#endif
