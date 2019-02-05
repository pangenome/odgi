#ifndef DANKGRAPH_HASH_MAP_HPP_INCLUDED
#define DANKGRAPH_HASH_MAP_HPP_INCLUDED

#include <cstdint>
#include <tuple>
#include <type_traits>
#include "spp.h"

namespace dg {

// Thomas Wang's integer hash function. In many implementations, std::hash
// is identity function for integers, which leads to performance issues.

inline size_t wang_hash_64(size_t key) {
    key = (~key) + (key << 21); // key = (key << 21) - key - 1;
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8); // key * 265
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4); // key * 21
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

// We need this second type for enable_if-based specialization
template<typename T, typename ImplementationMatched = void>
struct wang_hash;

// We can hash pointers
template<typename T>
struct wang_hash<T*> {
    size_t operator()(const T* pointer) const {
        return wang_hash_64(reinterpret_cast<size_t>(pointer));
    }
};

// We can hash any integer that can be implicitly widened to size_t.
// This covers 32 bit ints (which we need to be able to hash on Mac) and 64 bit ints
// This also coveres bools.
// See <https://stackoverflow.com/a/42679086>
template<typename T>
struct wang_hash<T, typename std::enable_if<std::is_integral<T>::value>::type> {
    size_t operator()(const T& x) const {
        static_assert(sizeof(T) <= sizeof(size_t), "widest hashable type is size_t");
        return wang_hash_64(static_cast<size_t>(x));
    }
};

// We can hash pairs
template<typename A, typename B>
struct wang_hash<std::pair<A, B>> {
    size_t operator()(const std::pair<A, B>& x) const {
        size_t hash_val = wang_hash<A>()(x.first);
        hash_val ^= wang_hash<B>()(x.second) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
        return hash_val;
    }
};


// Replacements for std::unordered_map.

template<typename K, typename V>
class hash_map : public spp::sparse_hash_map<K, V, wang_hash<K>> { };

template<typename K, typename V>
class string_hash_map : public spp::sparse_hash_map<K, V> { };

template<typename K, typename V>
class pair_hash_map : public spp::sparse_hash_map<K, V, wang_hash<K>> { };

template<typename K, typename V>
class hash_map<K*, V> : public spp::sparse_hash_map<K*, V, wang_hash<K*>> { };

// Replacements for std::unordered_set.

template<typename K>
class hash_set : public spp::sparse_hash_set<K, wang_hash<K>> { };

template<typename K>
class string_hash_set : public spp::sparse_hash_set<K> { };

template<typename K>
class pair_hash_set : public spp::sparse_hash_set<K, wang_hash<K>> { };

template<typename K>
class hash_set<K*> : public spp::sparse_hash_set<K*, wang_hash<K*>> { };

}

#endif
