//
//  dynamic_structs.cpp
//  
//

#include "dynamic_structs.hpp"


namespace dankgraph {
    
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

size_t SuccinctSplayTree::first_lower(const size_t& key) const {
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
        if (get_parent(y) != z) {
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
        set_right(y, x);
    }
    set_parent(x, y);
}
    
void SuccinctSplayTree::splay(size_t x) {
    while (get_parent(x) != 0) {
        if (get_parent(get_parent(x)) == 0) {
            if (get_left(get_parent(x)) == x) {
                right_rotate(get_parent(x));
            }
            else {
                left_rotate(get_parent(x));
            }
        }
        else if (get_left(get_parent(x)) == x && get_left(get_parent(get_parent(x))) == get_parent(x)) {
            right_rotate(get_parent(get_parent(x)));
            right_rotate(get_parent(x));
        }
        else if (get_right(get_parent(x)) == x && get_right(get_parent(get_parent(x))) == get_parent(x)) {
            left_rotate(get_parent(get_parent(x)));
            left_rotate(get_parent(x));
        }
        else if (get_left(get_parent(x)) == x && get_right(get_parent(get_parent(x))) == get_parent(x)) {
            right_rotate(get_parent(x));
            left_rotate(get_parent(x));
        }
        else {
            left_rotate(get_parent(x));
            right_rotate(get_parent(x));
        }
    }
}

void SuccinctSplayTree::replace(size_t u, size_t v ) {
    size_t u_parent = get_parent(u);
    if (u_parent == 0) {
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
    // return 1-based index
    return tree.size() / NODE_SIZE;
}
    
void SuccinctSplayTree::delete_node(size_t x) {
    
    // 1-based index of final node in the vector
    size_t last = tree.size() / NODE_SIZE;
    
    if (x != last) {
        
        // replace the values in x with the last node's values
        set_key(x, get_key(last));
        set_value(x, get_value(last));
        set_left(x, get_left(last));
        set_right(x, get_right(last));
        set_parent(x, get_parent(last));
        
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
    }
        
    // remove the last node from the vector
    for (size_t i = 0; i < NODE_SIZE; i++) {
        tree.pop();
    }
}
    
void SuccinctSplayTree::print_topology(std::ostream& out) const {
    std::function<void(std::ostream&,size_t, int, bool)> internal = [&](std::ostream& out, size_t x, int depth, bool right_child) {
        if (!right_child) {
            for (int i = 0; i < depth - 1; i++) {
                out << "\t\t";
            }
        }
        
        out << "\t->\t";
        
        if (x == 0) {
            out << "." << std::endl;
        }
        else {
            out << get_key(x) << " " << get_value(x);
            internal(out, get_right(x), depth + 1, true);
            internal(out, get_left(x), depth + 1, false);
        }
    };
    
    if (root == 0) {
        out << "." << std::endl;
    }
    else {
        out << get_key(root) << " " << get_value(root);
        internal(out, get_right(root), 1, true);
        internal(out, get_left(root), 1, false);
    }
}
    
void SuccinctSplayTree::print_vector(std::ostream& out) const {
    out << "root: " << root << std::endl;
    for (size_t i = 0; i < tree.size(); i++) {
        if (i % NODE_SIZE == 0 && i != 0) {
            out << "| ";
        }
        if (i % NODE_SIZE == 0) {
            out << "(" << i / NODE_SIZE + 1 << "): ";
        }
        out << tree.get(i) << " ";
    }
    out << std::endl;
}
}
