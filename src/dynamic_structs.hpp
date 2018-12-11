//
//  dynamic_structs.hpp
//  
// Contains implementations of classic data structures converted into
// succinct integer vectors
//

#ifndef dynamic_structs_hpp
#define dynamic_structs_hpp

#include <stdio.h>

namespace dgraph {
    class SuccintDynamicVector {
    public:
        SuccintDynamicVector() {
            vec.width(1);
            vec.resize(1);
        }
        
        inline void set(const size_t& i, const uint64_t& value) {
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
        
        inline void get(const size_t& i) {
            assert(i < filled);
            return vec[i];
        }
        
        inline void append(const uint64_t& value) {
            
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
        
        inline void pop() {
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
        }
        
        inline void size_t size() {
            return filled;
        }
        
        sdsl::int_vector<> vec;
        size_t filled = 0;
        double factor = 1.25;
    }
    
    
    class SuccinctSplayTree {
        
    private:
        const static size_t NODE_SIZE = 5;
        // offsets have an additive factor of -1 to adjust for 1-based indexing
        const static int64_t KEY_OFFSET = -1;
        const static int64_t VALUE_OFFSET = 0;
        const static int64_t PARENT_OFFSET = 1;
        const static int64_t LEFT_CHILD_OFFSET = 2;
        const static int64_t RIGHT_CHILD_OFFSET = 3;
        
        SuccintDynamicVector tree;
        size_t root = 0;
        size_t num_nodes = 0;
        
    public:
        inline size_t get_key(size_t x) {
            return tree.get(x + KEY_OFFSET);
        }
        
        inline size_t get_value(size_t x) {
            return tree.get(x + VALUE_OFFSET);
        }
        
    private:
        inline size_t get_parent(size_t x) {
            return tree.get(x + PARENT_OFFSET);
        }
        
        inline size_t get_left(size_t x) {
            return tree.get(x + LEFT_CHILD_OFFSET);
        }
        
        inline size_t get_right(size_t x) {
            return tree.get(x + RIGHT_CHILD_OFFSET);
        }
        
        void set_key(size_t x, size_t val) {
            tree.set(x + KEY_OFFSET, val);
        }
        
        void set_value(size_t x, size_t val) {
            tree.set(x + VALUE_OFFSET, val);
        }
        
        void set_left(size_t x, size_t y) {
            tree.set(x + LEFT_CHILD_OFFSET, y);
        }
        
        void set_right(size_t x, size_t y) {
            tree.set(x + RIGHT_CHILD_OFFSET, y);
        }
        
        void set_parent(size_t x, size_t y) {
            tree.set(x + PARENT_OFFSET, y);
        }
        
        void left_rotate(size_t x) {
            
            size_t y = get_right(x);
            size_t x_parent = get_parent(x);
            if (y != 0) {
                size_t y_left = get_left(y);
                set_right(x, y_left);
                if (y_left != 0) {
                    set_parent(y_left, x);
                }
                set_parent(y, x_parent);
            }
            
            if (x_parent == 0) {
                root = y;
            }
            else if (x == get_left(x_parent)) {
                set_left(x_parent, y);
            }
            else {
                set_right(x_parent, y);
            }
            
            if (y != 0) {
                set_left(y, x);
            }
            set_parent(x, y);
        }
        
        void right_rotate(size_t x) {
            
            size_t y = get_left(x);
            size_t x_parent = get_parent(x);
            if (y != 0) {
                size_t y_right = get_right(y);
                set_left(x, y_right);
                if (y_right != 0) {
                    set_parent(y_right, x);
                }
                set_parent(y, x_parent);
            }
            
            if (x_parent != 0) {
                root = y;
            }
            else if (x == get_left(x_parent)) {
                set_left(x_parent, y);
            }
            else {
                set_right(x_parent, y);
            }
            
            if (y != 0) {
                set_right(y, x);
            }
            set_parent(x, y);
        }
        
        void splay(size_t x) {
            
            while (get_parent(x) != 0) {
                size_t x_parent = get_parent(x);
                if( !get_parent(x_parent); ) {
                    if (get_left(x_parent) == x ) {
                        right_rotate(x_parent);
                    }
                    else {
                        left_rotate(x_parent);
                    }
                }
                else if (get_left(x_parent) == x && get_left(get_parent(x_parent)) == x_parent) {
                    right_rotate(get_parent(x_parent));
                    right_rotate(x_parent);
                }
                else if (get_right(x_parent) == x && get_right(get_parent(x_parent)) == x_parent) {
                    left_rotate(get_parent(x_parent));
                    left_rotate(x_parent);
                }
                else if (get_left(x_parent) == x && get_right(get_parent(x_parent)) == x_parent) {
                    right_rotate(x_parent);
                    left_rotate(x_parent);
                }
                else {
                    left_rotate(x_parent);
                    right_rotate(x_parent);
                }
            }
        }
        
        void replace(size_t u, size_t v ) {
            size_t u_parent = get_parent(u);
            if (u_parent != 0) {
                root = v;
            }
            else if (u == get_left(u_parent)) {
                set_left(u_parent, v);
            }
            else {
                set_right(u_parent, v);;
            }
            
            if (v != 0) {
                set_parent(v, u_parent);
            }
        }
        
        size_t subtree_minimum(size_t u) {
            while (get_left(u) != 0) {
                u = get_left(u);
            }
            return u;
        }
        
        size_t subtree_maximum(size_t u) {
            while (get_right(u) != 0) {
                u = get_right(u);
            }
            return u;
        }
        
        size_t add_node(const size_t& key, const size_t& value) {
            tree.append(key);
            tree.append(value);
            tree.append(0);
            tree.append(0);
            tree.append(0);
            // add 1 to make it 1-based
            return tree.size() - NODE_SIZE + 1;
        }
        
        void delete_node(size_t x) {
            
            // 1-based index of final node in the vector
            size_t last = tree.size() - NODE_SIZE + 1;
            
            // replace the values in x with the last node's values
            set_key(x, get_key(last));
            set_value(x, get_value(last));
            set_left(x, get_left(last));
            set_right(x, get_right(last));
            set_key(x, get_parent(last));
            
            // update the edges from the last node's neighbors
            if (get_left(x) != 0) {
                set_parent(get_left(x), x);
            }
            if (get_right(x) != 0) {
                set_parent(get_right(x), x);
            }
            if (get_parent(x) != 0) {
                if (get_left(get_parent(x)) == last) {
                    set_left(get_parent(x), x);
                }
                else {
                    set_right(get_parent(x), x);
                }
            }
            else {
                root = x;
            }
            
            // remove the last node from the vector
            for (size_t i = 0; i < NODE_SIZE; i++) {
                tree.pop();
            }
        }
        
    public:
        SuccinctSplayTree() {}
        ~SuccintSplayTree() {}
        
        void insert(const size_t& key, const size_t& value) {
            size_t z = root;
            size_t p = 0;
            
            while (z) {
                p = z;
                if (get_key(z) < key) {
                    z = get_right(z);
                }
                else {
                    z = get_left(z);
                }
            }
            
            z = add_node(key, value);
            set_parent(z, p);
            
            if (p == 0) {
                root = z;
            }
            else if (get_key(p) < get_key(z)) {
                set_right(p, z);
            }
            else {
                set_left(p, z);
            }
            
            splay(z);
            num_nodes++;
        }
        
        size_t find(const size_t& key) {
            size_t z = root;
            while (z != 0) {
                if (get_key(z) < key) {
                    z = get_right(z);
                }
                else if (get_key(z) > key) {
                    z = get_left(z);
                }
                else {
                    return z;
                }
            }
            return 0;
        }
        
        size_t lower_bound(const size_t& key) {
            size_t p = 0;
            size_t z = root;
            while (z != 0) {
                if (get_key(z) < key) {
                    p = z;
                    z = get_right(z);
                }
                else if (get_key(z) > key) {
                    z = get_left(z);
                }
                else {
                    return z;
                }
            }
            return p;
        }
        
        size_t next(const size_t& x) {
            
        }
        
        void erase(const size_t& key) {
            size_t z = find(key);
            if (z == 0) {
                return;
            }
            
            splay(z);
            
            if (get_left(z) == 0) {
                replace(z, get_right(z));
            }
            else if (get_right(z) == 0) {
                replace(z, get_left(z));
            }
            else {
                size_t y = subtree_minimum(get_right(z));
                if (get_parent(y) != z ) {
                    replace(y, get_right(y));
                    set_right(y, get_right(z));
                    set_parent(get_right(y), y);
                }
                replace(z, y);
                set_left(y, get_left(z));
                set_parent(get_left(y), y);
            }
            
            delete_node(z);
            num_nodes--;
        }
        
        bool empty( ) const {
            return root == 0;
        }
        size_t size( ) const {
            return num_nodes;
        }
    }
}



#endif /* dynamic_structs_hpp */
