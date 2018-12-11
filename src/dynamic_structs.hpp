//
//  dynamic_structs.hpp
//  
// Contains implementations of classic data structures converted into
// succinct integer vectors
//

#ifndef dynamic_structs_hpp
#define dynamic_structs_hpp

#include <stdio.h>

namespace dankgraph {
    class SuccinctDynamicVector {
    public:
        SuccinctDynamicVector();
        ~SuccinctDynamicVector();
        
        inline void set(const size_t& i, const uint64_t& value);
        
        inline void get(const size_t& i) const;
        
        inline void append(const uint64_t& value);
        
        inline void pop();
        
        inline size_t size() const;
        
        inline bool empty() const;
        
    private:
        
        sdsl::int_vector<> vec;
        size_t filled = 0;
        double factor = 1.25;
    }
    
    
    class SuccinctSplayTree {
        
    public:
        SuccinctSplayTree();
        ~SuccintSplayTree();
        
        void insert(const size_t& key, const size_t& value);
        
        void erase(const size_t& key);
        
        bool empty() const;
        
        size_t size() const;
        
        size_t find(const size_t& key);
        
        size_t lower_bound(const size_t& key);
        
        size_t next(const size_t& x);
        
        inline size_t get_key(size_t x);
        
        inline size_t get_value(size_t x);
        
    private:
        const static size_t NODE_SIZE = 5;
        // offsets have an additive factor of -1 to adjust for 1-based indexing
        const static int64_t KEY_OFFSET = -1;
        const static int64_t VALUE_OFFSET = 0;
        const static int64_t PARENT_OFFSET = 1;
        const static int64_t LEFT_CHILD_OFFSET = 2;
        const static int64_t RIGHT_CHILD_OFFSET = 3;
        
        SuccinctDynamicVector tree;
        size_t root = 0;
        size_t num_nodes = 0;
        
        inline size_t get_parent(size_t x);
        
        inline size_t get_left(size_t x);
        
        inline size_t get_right(size_t x);
        
        inline void set_key(size_t x, size_t val);
        
        inline void set_value(size_t x, size_t val);
        
        inline void set_left(size_t x, size_t y);
        
        inline void set_right(size_t x, size_t y);
        
        inline void set_parent(size_t x, size_t y);
        
        void left_rotate(size_t x);
        
        void right_rotate(size_t x);
        
        void splay(size_t x);
        
        void replace(size_t u, size_t v );
        
        size_t subtree_minimum(size_t u);
        
        size_t subtree_maximum(size_t u);
        
        size_t add_node(const size_t& key, const size_t& value);
        
        void delete_node(size_t x);
    }
    
    
    
    
    /// Inline functions
    
    inline void SuccinctDynamicVector::set(const size_t& i, const uint64_t& value) {
        assert(i < filled);
        
        uint8_t width = vec.width();
        uint64_t mask = numeric_limits<uint64_t>::max() << width;
        while (mask & value) {
            width++;
            mask = numeric_limits<uint64_t>::max() << width;
        }
        
        if (width > vec.width()) {
            vec.width(width)
        }
        
        vec[i] = value;
    }
    
    inline void SuccinctDynamicVector::get(const size_t& i) const {
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
            vec = move(tmp);
        }
        
        filled++;
        set(filled - 1, value);
    }
    
    inline void SuccinctDynamicVector::pop() {
        filled--;
        size_t shrink_capacity = vec.size() / (factor * factor)
        if (filled <= shrink_capacity) {
            sdsl::int_vector<> tmp;
            tmp.width(vec.width());
            tmp.resize(shrink_capacity);
            for (size_t i = 0; i < vec.size(); i++) {
                tmp[i] = vec[i];
            }
            vec = move(tmp);
        }
    }'
    
    inline size_t SuccinctDynamicVector::size() const {
        return filled;
    }'
    
    inline bool SuccinctDynamicVector::empty() const {
        return filled == 0;
    }
    
    inline bool SuccinctSplayTree::empty( ) const {
        return root == 0;
    }
    
    inline size_t SuccinctSplayTree::size( ) const {
        return num_nodes;
    }
    
    inline size_t SuccinctSplayTree::get_key(size_t x) {
        return tree.get(x + KEY_OFFSET);
    }
    
    inline size_t SuccinctSplayTree::get_value(size_t x) {
        return tree.get(x + VALUE_OFFSET);
    }
    
    inline size_t SuccinctSplayTree::get_parent(size_t x) {
        return tree.get(x + PARENT_OFFSET);
    }
    
    inline size_t SuccinctSplayTree::get_left(size_t x) {
        return tree.get(x + LEFT_CHILD_OFFSET);
    }
    
    inline size_t SuccinctSplayTree::get_right(size_t x) {
        return tree.get(x + RIGHT_CHILD_OFFSET);
    }
    
    inline void SuccinctSplayTree::set_key(size_t x, size_t val) {
        tree.set(x + KEY_OFFSET, val);
    }
    
    inline void SuccinctSplayTree::set_value(size_t x, size_t val) {
        tree.set(x + VALUE_OFFSET, val);
    }
    
    inline void SuccinctSplayTree::set_left(size_t x, size_t y) {
        tree.set(x + LEFT_CHILD_OFFSET, y);
    }
    
    inline void SuccinctSplayTree::set_right(size_t x, size_t y) {
        tree.set(x + RIGHT_CHILD_OFFSET, y);
    }
    
    inline void SuccinctSplayTree::set_parent(size_t x, size_t y) {
        tree.set(x + PARENT_OFFSET, y);
    }
}



#endif /* dynamic_structs_hpp */
