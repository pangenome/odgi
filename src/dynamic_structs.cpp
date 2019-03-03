//
//  dynamic_structs.cpp
//  
//

#include "dynamic_structs.hpp"


namespace odgi {
    
const double SuccinctDynamicVector::factor = 1.25;
    
SuccinctDynamicVector::SuccinctDynamicVector() {
    vec.width(1); // by default we start as a bitvector
    vec.resize(1);
}
    
SuccinctDynamicVector::~SuccinctDynamicVector() {
    
}
    
SuccinctSplayTree::SuccinctSplayTree() {
    
}
    
SuccinctSplayTree::~SuccinctSplayTree() {
    
}
    
void SuccinctSplayTree::insert(const size_t& key, const size_t& value) {
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

size_t SuccinctSplayTree::find(const size_t& key) const {
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

size_t SuccinctSplayTree::lower_bound(const size_t& key) const {
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

size_t SuccinctSplayTree::next(const size_t& x) const {
    if (get_right(x) != 0) {
        return subtree_minimum(get_right(x));
    }
    else {
        size_t c = x;
        size_t z = get_parent(x);
        while (z != 0) {
            if (c == get_left(z)) {
                break;
            }
                
            c = z;
            z = get_parent(z);
        }
        return z;
    }
}
    
void SuccinctSplayTree::erase(const size_t& key) {
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
    
void SuccinctSplayTree::left_rotate(size_t x) {
        
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
    
void SuccinctSplayTree::right_rotate(size_t x) {
        
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
    
void SuccinctSplayTree::splay(size_t x) {
    while (get_parent(x) != 0) {
        size_t x_parent = get_parent(x);
        if (get_parent(x_parent) == 0) {
            if (get_left(x_parent) == x) {
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

void SuccinctSplayTree::replace(size_t u, size_t v ) {
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

size_t SuccinctSplayTree::subtree_minimum(size_t u) const {
    while (get_left(u) != 0) {
        u = get_left(u);
    }
    return u;
}
    
size_t SuccinctSplayTree::subtree_maximum(size_t u) const {
    while (get_right(u) != 0) {
        u = get_right(u);
    }
    return u;
}
    
size_t SuccinctSplayTree::add_node(const size_t& key, const size_t& value) {
    tree.append(key);
    tree.append(value);
    tree.append(0);
    tree.append(0);
    tree.append(0);
    // add 1 to make it 1-based
    return tree.size() - NODE_SIZE + 1;
}
    
void SuccinctSplayTree::delete_node(size_t x) {
        
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
}
