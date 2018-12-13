//
//  dynamic_structs.hpp
//  
// Contains implementations of classic data structures converted into
// succinct integer vectors
//

#ifndef dynamic_structs_hpp
#define dynamic_structs_hpp

#include <cstdio>
#include <cstdint>
#include "sdsl/bit_vectors.hpp"
#include "sdsl/enc_vector.hpp"
#include "sdsl/dac_vector.hpp"
#include "sdsl/vlc_vector.hpp"
#include "sdsl/wavelet_trees.hpp"
#include "sdsl/csa_wt.hpp"
#include "sdsl/suffix_arrays.hpp"

namespace dankgraph {
    
/*
 * A dynamic integer vector that maintains integers in bit-compressed form.
 * Automatically adjusts bit-width for entries depending on input data.
 */
class SuccinctDynamicVector {
public:
    /// Constructor (starts empty)
    SuccinctDynamicVector();
        
    /// Destructor
    ~SuccinctDynamicVector();
        
    /// Set the i-th value
    inline void set(const size_t& i, const uint64_t& value);
        
    /// Returns the i-th value
    inline uint64_t get(const size_t& i) const;
        
    /// Add a value to the end
    inline void append(const uint64_t& value);
        
    /// Remove the last value
    inline void pop();
        
    /// Returns the number of values
    inline size_t size() const;
        
    /// Returns true if there are no entries and false otherwise
    inline bool empty() const;

    /// Clears the backing vector
    inline void clear();
        
private:
        
    // the underlying vector representation
    sdsl::int_vector<> vec;
    // tracker for number of values
    size_t filled = 0;
    // geometric expansion factor
    static const double factor;
};
    
    
class SuccinctSplayTree {

public:
    SuccinctSplayTree(void);
    ~SuccinctSplayTree(void);
    
    /// Insert a key-value pair. If the key already exists, the current value
    /// will be replaced with the given value.
    void insert(const size_t& key, const size_t& value);
    
    /// Erase the key-value pair associated with the key. If the key does not
    /// exist, do nothing.
    void erase(const size_t& key);
    
    /// Returns true if there are no entries, otherwise false.
    bool empty() const;
    
    /// Returns the number of entries.
    size_t size() const;
    
    /// Returns a handle to the key-value pair associated with a key, or 0 if
    /// there the key does not exist.
    size_t find(const size_t& key) const;
    
    /// Returns a handle to the key-value pair with the largest key that is less-than
    /// or equal to the given key, or 0 if the given key is less than the minimum or
    /// the tree is empty.
    size_t first_lower(const size_t& key) const;
    
    /// Returns the handle to the key-value pair with the next-smallest key to the
    /// given handle.
    size_t next(const size_t& x) const;
    
    /// Returns the key of a handle.
    inline size_t get_key(const size_t& x) const;
    
    /// Returns the value of a handle.
    inline size_t get_value(const size_t& x) const;

        
private:
    const static size_t NODE_SIZE = 5;
    const static int64_t KEY_OFFSET = 0;
    const static int64_t VALUE_OFFSET = 1;
    const static int64_t PARENT_OFFSET = 2;
    const static int64_t LEFT_CHILD_OFFSET = 3;
    const static int64_t RIGHT_CHILD_OFFSET = 4;
        
    SuccinctDynamicVector tree;
    size_t root = 0;
    size_t num_nodes = 0;
        
    inline size_t get_parent(size_t x) const;
    
    inline size_t get_left(size_t x) const;
    
    inline size_t get_right(size_t x) const;
        
    inline void set_key(size_t x, size_t val);
        
    inline void set_value(size_t x, size_t val);
        
    inline void set_left(size_t x, size_t y);
        
    inline void set_right(size_t x, size_t y);
        
    inline void set_parent(size_t x, size_t y);
        
    void left_rotate(size_t x);
        
    void right_rotate(size_t x);
        
    void splay(size_t x);
        
    void replace(size_t u, size_t v );
        
    size_t subtree_minimum(size_t u) const;
    
    size_t subtree_maximum(size_t u) const;
    
    size_t add_node(const size_t& key, const size_t& value);
        
    void delete_node(size_t x);
    
    void print_topology(std::ostream& out) const;
    void print_vector(std::ostream& out) const;
};
    
    
    
    
/// Inline functions
    
inline void SuccinctDynamicVector::set(const size_t& i, const uint64_t& value) {
    assert(i < filled);
        
    uint8_t width = vec.width();
    uint64_t mask = std::numeric_limits<uint64_t>::max() << width;
    while (mask & value) {
        width++;
        mask = std::numeric_limits<uint64_t>::max() << width;
    }
        
    if (width > vec.width()) {
        sdsl::int_vector<> wider_vec;
        wider_vec.width(width);
        wider_vec.resize(vec.size());
        for (size_t i = 0; i < filled; i++) {
            wider_vec[i] = vec[i];
        }
        vec = std::move(wider_vec);
    }
        
    vec[i] = value;
}
    
inline uint64_t SuccinctDynamicVector::get(const size_t& i) const {
    assert(i < filled);
    return vec[i];
}
    
inline void SuccinctDynamicVector::append(const uint64_t& value) {
        
    if (filled == vec.size()) {
        size_t new_capacity = size_t(vec.size() * factor) + 1;
        sdsl::int_vector<> tmp;
        tmp.width(vec.width());
        tmp.resize(new_capacity);
        for (size_t i = 0; i < filled; i++) {
            tmp[i] = vec[i];
        }
        vec = std::move(tmp);
    }
        
    filled++;
    set(filled - 1, value);
}
    
inline void SuccinctDynamicVector::pop() {
    filled--;
    size_t shrink_capacity = vec.size() / (factor * factor);
    if (filled < shrink_capacity) {
        sdsl::int_vector<> tmp;
        tmp.width(vec.width());
        tmp.resize(shrink_capacity);
        for (size_t i = 0; i < filled; i++) {
            tmp[i] = vec[i];
        }
        vec = std::move(tmp);
    }
}
    
inline size_t SuccinctDynamicVector::size() const {
    return filled;
}
    
inline bool SuccinctDynamicVector::empty() const {
    return filled == 0;
}

inline void SuccinctDynamicVector::clear() {
    vec.resize(0);
    filled = 0;
}
    
inline bool SuccinctSplayTree::empty() const {
    return root == 0;
}
    
inline size_t SuccinctSplayTree::size() const {
    return num_nodes;
}
    
inline size_t SuccinctSplayTree::get_key(const size_t& x) const {
    return tree.get((x - 1) * NODE_SIZE + KEY_OFFSET);
}

inline size_t SuccinctSplayTree::get_value(const size_t& x) const {
    return tree.get((x - 1) * NODE_SIZE + VALUE_OFFSET);
}

inline size_t SuccinctSplayTree::get_parent(size_t x) const {
    return tree.get((x - 1) * NODE_SIZE + PARENT_OFFSET);
}

inline size_t SuccinctSplayTree::get_left(size_t x) const {
    return tree.get((x - 1) * NODE_SIZE + LEFT_CHILD_OFFSET);
}

inline size_t SuccinctSplayTree::get_right(size_t x) const {
    return tree.get((x - 1) * NODE_SIZE + RIGHT_CHILD_OFFSET);
}
    
inline void SuccinctSplayTree::set_key(size_t x, size_t val) {
    tree.set((x - 1) * NODE_SIZE + KEY_OFFSET, val);
}
    
inline void SuccinctSplayTree::set_value(size_t x, size_t val) {
    tree.set((x - 1) * NODE_SIZE + VALUE_OFFSET, val);
}
    
inline void SuccinctSplayTree::set_left(size_t x, size_t y) {
    tree.set((x - 1) * NODE_SIZE + LEFT_CHILD_OFFSET, y);
}
    
inline void SuccinctSplayTree::set_right(size_t x, size_t y) {
    tree.set((x - 1) * NODE_SIZE + RIGHT_CHILD_OFFSET, y);
}
    
inline void SuccinctSplayTree::set_parent(size_t x, size_t y) {
    tree.set((x - 1) * NODE_SIZE + PARENT_OFFSET, y);
}
}



#endif /* dynamic_structs_hpp */
