#include "kmer.hpp"

namespace odgi {

namespace kmer {

uint64_t seq2bit(const char* s, uint64_t k, bool& is_dna) {
    is_dna = true;
    uint64_t v = 0;
    for (uint64_t i = 0; i < k; ++i) {
        v <<= 2;
        switch (*(s+i)) {
        case 'A': case 'a':
            // 0b00, default case
            break;
        case 'T': case 't':
            v |= 0b01;
            break;
        case 'G': case 'g':
            v |= 0b10;
            break;
        case 'C': case 'c':
            v |= 0b11;
            break;
        default:
            is_dna = false;
            break;
        }
    }
    return v;
}

std::string unseq2bit(uint64_t v, const uint64_t& size) {
    assert(size <= sizeof(uint64_t)*4);
    std::string s(size, '\0');
    for (int64_t i = size-1; i >= 0; --i) {
        switch (v & 0b11) {
        case 0b00:
            s[i] = 'A';
            break;
        case 0b01:
            s[i] = 'T';
            break;
        case 0b10:
            s[i] = 'G';
            break;
        case 0b11:
            s[i] = 'C';
            break;
        default:
            break;
        }
        v >>= 2;
    }
    return s;
}

double entropy(uint64_t v, const uint64_t& size) {
    uint64_t counts[4] = {0, 0, 0, 0};
    for (int64_t i = size-1; i >= 0; --i) {
        ++counts[v & 0b11];
        v >>= 2;
    }
    double ent = 0;
    double ln2 = log(2);
    for (uint8_t i = 0; i < 4; ++i) {
        if (counts[i]) {
            double f = (double)counts[i]/size;
            ent += f * log(f)/ln2;
        }
    }
    return -ent;
}

double gc_rate(uint64_t v, const uint64_t& size) {
    uint64_t gc = 0;
    for (int64_t i = size-1; i >= 0; --i) {
        switch (v & 0b11) {
        case 0b10: // G
        case 0b11: // C
            ++gc;
            break;
        default:
            break;
        }
        v >>= 2;
    }
    return (double)gc/size;
}

}

}
