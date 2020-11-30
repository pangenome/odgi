#pragma once

#include <iostream>
#include <cstdint>
#include <string>
#include <handlegraph/util.hpp>
#include <vector>
#include <map>
#include <utility>
#include <cstring>
#include <cassert>
#include "bmap.hpp"
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
    bmap::bmap<uint64_t, uint64_t> encoding;
    dyn::hacked_vector decoding;
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

    uint64_t encode(const uint64_t& id);
    uint64_t decode(const uint64_t& id) const;

    uint64_t sequence_size(void) const;
    const std::string& get_sequence(void) const;
    void set_sequence(const std::string& seq);
    void for_each_edge(const std::function<bool(nid_t other_id,
                                                bool on_rev,
                                                bool other_rev,
                                                bool to_curr)>& func) const;
    void add_edge(const uint64_t& other_id,
                  const bool& on_rev,
                  const bool& to_rev,
                  const bool& to_curr);
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
    void remove_path_step(const uint64_t& rank);
    void update_path_last_bytes(void);
    void clear(void);
    void clear_edges(void);
    void clear_paths(void);
    void clear_encoding(void);
    uint64_t serialize(std::ostream& out) const;
    void load(std::istream& in);
    void display(void) const;

};

}
