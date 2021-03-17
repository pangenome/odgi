//----------------------------------------------------------------------
// Authors: Michael Calmette, Erik Garrison
// File:    bmap.hpp
//
// A map data stucture based on a sorted vector of key/value pairs.
// binary search with an add, remove, find and sort functions.
//
//----------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <vector>
#include "dynamic.hpp"

namespace bmap {

template <typename K, typename V> class bmap {

  public:
    // add a new key-value pair into the collection
    void add(const K &a_key, const V &a_val);

    // remove a key-value pair from the collectiona
    void remove(const K &a_key);

    // find and return the value associated with the key
    bool find(const K &search_key, V &the_val) const;

    // find and return the values with keys >= to k1 and <= to k2
    void find(const K &k1, const K &k2, std::vector<V> &vals) const;

    // return all of the keys in the collection
    //void keys(std::vector<K> &all_keys) const;

    // return all of the keys in ascending (sorted) order
    //void sort(std::vector<K> &all_keys_sorted) const;

    // return the number of key-value pairs in the collection
    int size() const;

    // save to the stream
    uint64_t serialize(std::ostream& out) const;

    // load from the stream
    void load(std::istream& in);

    // clear storage
    void clear();

private:
    // helper function for binary search
    bool binsearch(const K &key, int &index) const {
        if (kv_list.size() == 0)
            return false;

        int s;
        // markers to cut list in half
        int end = kv_list.size()/2 - 1;
        int begin = 0;

        // check if a key is there, find index pos
        while (end >= begin) {
            s = begin + ((end - begin) / 2);
            if (kv_list.at(s*2) == key) {
                index = s;
                return true;
            } else if (kv_list.at(s*2) < key)
                begin = s + 1;

            else
                end = s - 1;
        }

        // key is not there, find index pos
        //s = 0;
        for (s = 0; s < kv_list.size()/2; ) { //std::pair<K, V> p : kv_list) {
            if (key < kv_list.at(s*2)) {
                index = s;
                return false;
            }
            ++s;
        }
        index = s; // end of list
        return false;
    };

    // vector storage
    //std::vector<std::pair<K, V>> kv_list;
    dyn::hacked_vector kv_list;
};

// TODO: Implement functions here ...

// add key value pair to vector
template <typename K, typename V>
void bmap<K, V>::add(const K &a_key, const V &a_val) {
    std::pair<K, V> p(a_key, a_val);
    int idx = 0;
    // get index position to add key
    bool check = binsearch(a_key, idx);
    if (size() == 0) {
        kv_list.push_back(p.first); // empty vector
        kv_list.push_back(p.second);
        return;
    }
    kv_list.insert(idx*2, p.second);
    kv_list.insert(idx*2, p.first);
}

// remove key from vector
template <typename K, typename V> void bmap<K, V>::remove(const K &a_key) {
    int idx = 0;
    // check for key in vector
    bool check = binsearch(a_key, idx);
    if (check == true) {
        kv_list.remove(idx*2);
        kv_list.remove(idx*2);
    }
}

// find key in vector
template <typename K, typename V>
bool bmap<K, V>::find(const K &search_key, V &the_val) const {
    int idx = 0;
    // check if in vector and index if it is
    bool check = binsearch(search_key, idx);
    if (check == true) {
        the_val = kv_list.at(idx*2+1); // set val
        return true;
    }
    return false;
}

// find and return vals >= k1 and <= k2
template <typename K, typename V>
void bmap<K, V>::find(const K &k1, const K &k2, std::vector<V> &vals) const {
    int idx1 = 0;
    int idx2 = 0;
    bool checkOne = binsearch(k1, idx1);
    bool checkTwo = binsearch(k2, idx2);

    if (checkOne == true && checkTwo == true) {
        for (int i = idx1; i <= idx2; i++) {
            vals.push_back(kv_list.at(i*2+1));
        }
    }
}

// return all keys in sorted order
/*
template <typename K, typename V>
void bmap<K, V>::keys(std::vector<K> &all_keys) const {
    for (std::pair<K, V> p : kv_list) {
        all_keys.push_back(p.first); // push keys
    }
}
*/

// sort keys
/*
template <typename K, typename V>
void bmap<K, V>::sort(std::vector<K> &all_keys_sorted) const {
    if (kv_list.size() == 0)
        return;
    for (std::pair<K, V> p : kv_list) {
        all_keys_sorted.push_back(p.first);
    }
    std::sort(all_keys_sorted.begin(), all_keys_sorted.end());
}
*/

// return vector size
template <typename K, typename V> int bmap<K, V>::size() const {
    return kv_list.size() / 2;
}


template <typename K, typename V>
uint64_t bmap<K, V>::serialize(std::ostream &out) const {
    return kv_list.serialize(out);
    /*
    uint64_t w_size = sizeof(std::pair<K,V>)*kv_list.size();
    if (w_size) {
        out.write(reinterpret_cast<const char*>(&w_size), sizeof(w_size));
        out.write(reinterpret_cast<const char*>(kv_list.data()), w_size);
    }
    return w_size;
    */
}

template <typename K, typename V>
void bmap<K, V>::load(std::istream &in) {
    kv_list.load(in);
    /*
    uint64_t w_size;
    in.read(reinterpret_cast<char*>(&w_size), sizeof(w_size));
    kv_list = std::vector<std::pair<K,V>>(w_size);
    if (w_size) {
        in.read(reinterpret_cast<char*>(kv_list.data()), w_size);
    }
    */
}

template <typename K, typename V>
void bmap<K, V>::clear() {
    //kv_list = std::vector<std::pair<K,V>>();
    dyn::hacked_vector null_iv;
    kv_list = null_iv;
}

} // namespace bmap
