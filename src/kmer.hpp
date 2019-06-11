#pragma once

#include <string>
#include <cassert>
#include <cmath>

namespace odgi {

namespace kmer {

uint64_t seq2bit(const char* s, uint64_t k, bool& is_dna);
std::string unseq2bit(uint64_t v, const uint64_t& size);
double entropy(uint64_t v, const uint64_t& size);
double gc_rate(uint64_t v, const uint64_t& size);

}

}
