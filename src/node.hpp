#pragma once

#include <iostream>
#include <cstdint>
#include <string>
#include <handlegraph/util.hpp>
#include <vector>
#include <map>
#include <utility>
#include <cstring>
#include "dynamic.hpp"
#include "varint.hpp"

namespace odgi {

using namespace handlegraph;
using nid_t = handlegraph::nid_t;
const uint8_t EDGE_RECORD_LENGTH = 2;
const uint8_t PATH_RECORD_LENGTH = 6;

/// A node object with the sequence, its edge lists, and paths
class node_t {
    std::string sequence;
    dyn::hacked_vector encoding;
    dyn::hacked_vector edges;
    dyn::hacked_vector paths;
public:
    /// edge type conversion
    /// 1 = fwd->fwd, 2 = fwd->rev, 3 = rev->fwd, 4 = rev->rev
    struct edge_helper {
        inline static uint8_t pack(bool other_rev, bool to_curr, bool on_rev) {
            return other_rev | (on_rev << 1) | (to_curr << 2);
        }
        inline static uint8_t unpack_on_rev(uint8_t edge) {
            return edge & (1 << 1);
        }
        inline static uint8_t unpack_other_rev(uint8_t edge) {
            return edge & 1;
        }
        inline static uint8_t unpack_to_curr(uint8_t edge) {
            return edge & (1 << 2);
        }
    };
    struct step_type_helper {
        inline static uint8_t pack(const bool& is_rev, const bool& is_start, const bool& is_end) {
            return is_rev | (is_start << 1) | (is_end << 2);
        }
        inline static bool unpack_is_rev(const uint8_t& type) {
            return type & 1;
        }
        inline static bool unpack_is_start(const uint8_t& type) {
            return type & (1 << 1);
        }
        inline static bool unpack_is_end(const uint8_t& type) {
            return type & (1 << 2);
        }
    };
    inline const uint64_t edge_count(void) const { return edges.size()/EDGE_RECORD_LENGTH; }
    inline const uint64_t path_count(void) const { return paths.size()/PATH_RECORD_LENGTH; }
    struct step_t {
        uint64_t path_id;
        bool is_rev;
        bool is_start;
        bool is_end;
        uint64_t prev_id;
        uint64_t prev_rank;
        uint64_t next_id;
        uint64_t next_rank;
    };
    /*
    struct step_t {
        uint64_t data[5] = { 0, 0, 0, 0, 0 }; // PATH_RECORD_LENGTH
        step_t(void) { }
        step_t(const uint64_t& id,
               const bool& rev,
               const bool& is_start,
               const bool& is_end,
               const uint64_t& prev_id,
               const uint64_t& prev_rank,
               const uint64_t& next_id,
               const uint64_t& next_rank) {
            data[0] = pack_step(id, rev, is_start, is_end);
            data[1] = prev_id;
            data[2] = prev_rank;
            data[3] = next_id;
            data[4] = next_rank;
        }
        inline const uint64_t path_id(void) const { return get_step_path_id(data[0]); }
        inline const bool is_rev(void) const { return get_step_is_rev(data[0]); }
        inline const bool is_start(void) const { return get_step_is_start(data[0]); }
        inline const bool is_end(void) const { return get_step_is_end(data[0]); }
        inline const uint64_t prev_id(void) const { return data[1]; }
        inline const uint64_t prev_rank(void) const { return data[2]; }
        inline const uint64_t next_id(void) const { return data[3]; }
        inline const uint64_t next_rank(void) const { return data[4]; }
        inline void set_path_step(const uint64_t& i) { data[0] = i; }
        inline void set_path_id(const uint64_t& i) { data[0] = pack_step(i, is_rev(), is_start(), is_end()); }
        inline void set_is_rev(const bool& b) { data[0] &= b << 2; }
        inline void set_is_start(const bool& b) { data[0] &= b << 1; }
        inline void set_is_end(const bool& b) { data[0] &= b; }
        inline void set_prev_id(const uint64_t& i) { data[1] = i; }
        inline void set_prev_rank(const uint64_t& i) { data[2] = i; }
        inline void set_next_id(const uint64_t& i) { data[3] = i; }
        inline void set_next_rank(const uint64_t& i) { data[4] = i; }
    };
    */
    // path step helpers
    /*
    inline static uint64_t get_step_path_id(const uint64_t& step) {
        return step >> 3;
    }
    inline static bool get_step_is_rev(const uint64_t& step) {
        return step & (1 << 2);
    }
    inline static bool get_step_is_start(const uint64_t& step) {
        return step & (1 << 1);
    }
    inline static bool get_step_is_end(const uint64_t& step) {
        return step & 1;
    }
    */

    uint64_t sequence_size(void) const;
    const std::string& get_sequence(void) const;
    //dyn::hacked_vector& get_edges(void) const;
    void set_sequence(const std::string& seq);
    //dyn::hacked_vector get_edges(void) const;
    void for_each_edge(const std::function<bool(nid_t other_id,
                                                bool on_rev,
                                                bool other_rev,
                                                bool to_curr)>& func) const;
    void add_edge(const uint64_t& other_id,
                  const bool& on_rev,
                  const bool& to_rev,
                  const bool& to_curr);
    //void remove_edge(const uint64_t& rank);
    bool remove_edge(const uint64_t& target_id,
                     const bool& target_rev,
                     const bool& ends_here,
                     const bool& is_rev);
    void add_path_step(const uint64_t& path_id, const bool& is_rev,
                       const bool& is_start, const bool& is_end,
                       const uint64_t& prev_id, const uint64_t& prev_rank,
                       const uint64_t& next_id, const uint64_t& next_rank);
    void add_path_step(const node_t::step_t& step);
    const step_t get_path_step(const uint64_t& rank) const;
    //const node_t::step_t node_t::get_path_step(const uint64_t& rank) const;
    const std::vector<step_t> get_path_steps(void) const;
    void set_path_step(const uint64_t& rank, const uint64_t& path_id, const bool& is_rev,
                       const bool& is_start, const bool& is_end,
                       const uint64_t& prev_id, const uint64_t& prev_rank,
                       const uint64_t& next_id, const uint64_t& next_rank);
    void set_path_step(const uint64_t& rank, const step_t& step);
    bool step_is_start(const uint64_t& rank) const;
    bool step_is_end(const uint64_t& rank) const;
    void for_each_path_step(const std::function<bool(step_t step)>& func) const;
    std::pair<std::map<uint64_t, std::pair<uint64_t, bool>>, // path fronts and backs
          std::map<uint64_t, std::pair<uint64_t, bool>>> flip_paths(void);
    //const step_t get_path_step(const uint64_t& rank) const;
    void remove_path_step(const uint64_t& rank);
    void update_path_last_bytes(void);
    void clear(void);
    void clear_edges(void);
    void clear_paths(void);
    uint64_t serialize(std::ostream& out) const;
    void load(std::istream& in);
    void display(void) const;

};

}
