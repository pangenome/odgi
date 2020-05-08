#pragma once

#include <cstdint>

namespace odgi {

uint32_t djb2_hash32(const char *str);
uint64_t djb2_hash64(const char *str);

}
