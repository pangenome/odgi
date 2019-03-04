#ifndef dgraph_node_hpp
#define dgraph_node_hpp

#include <iostream>
#include <cstdint>
#include <string>
#include <handlegraph/util.hpp>
#include "varint.hpp"

namespace odgi {

using namespace handlegraph;
using id_t = handlegraph::id_t;
const uint8_t EDGE_RECORD_LENGTH = 2;
const uint8_t PATH_RECORD_LENGTH = 5;

/// A node object with the sequence, its edge lists, and paths
class node_t {
    std::vector<uint8_t> bytes;
public:
    node_t(void) { };
    ~node_t(void) { };
    node_t(const node_t& other) {
        this->bytes = other.bytes;
    }
    node_t(node_t&& other) noexcept {
        this->bytes = other.bytes;
    }
    node_t& operator=(const node_t& other) {
        this->bytes = other.bytes;
        return *this;
    }
    node_t(const std::string& sequence) {
        init(sequence);
    }
    uint64_t size(void) {
        return bytes.size();
    }
    void clear(void) {
        bytes.clear();
    }
    uint8_t* data(void) {
        return bytes.data();
    }
    void resize(uint64_t i) {
        bytes.resize(i);
    }
    struct layout_t {
        uint64_t data[5];
        inline const uint64_t layout_bytes(void) const { return data[0]; }
        inline const uint64_t seq_start(void) const { return data[0]; }
        inline const uint64_t seq_bytes(void) const { return data[1]; }
        inline const uint64_t edge_start(void) const { return layout_bytes()+seq_bytes(); }
        inline const uint64_t edge_count(void) const { return data[2]; }
        inline const uint64_t edge_bytes(void) const { return data[3]; }
        inline const uint64_t path_start(void) const { return layout_bytes()+seq_bytes()+edge_bytes(); }
        inline const uint64_t path_count(void) const { return data[4]; }
        inline void set_layout_bytes(const uint64_t& i) { data[0] = i; }
        inline void set_seq_bytes(const uint64_t& i) { data[1] = i; }
        inline void set_edge_count(const uint64_t& i) { data[2] = i; }
        inline void set_edge_bytes(const uint64_t& i) { data[3] = i; }
        inline void set_path_count(const uint64_t& i) { data[4] = i; }
    };
    struct step_t {
        uint64_t data[5]; // PATH_RECORD_LENGTH
        inline const uint64_t path_id(void) const { return step_path_id(data[0]); }
        inline const bool is_rev(void) const { return step_is_rev(data[0]); }
        inline const uint64_t prev_id(void) const { return data[1]; }
        inline const uint64_t prev_rank(void) const { return data[2]; }
        inline const uint64_t next_id(void) const { return data[3]; }
        inline const uint64_t next_rank(void) const { return data[4]; }
        inline void set_path_step(const uint64_t& i) { data[0] = i; }
        inline void set_path_id(const uint64_t& i) { data[0] = pack_step(i, is_rev()); }
        inline void set_is_rev(const bool& b) { data[0] = pack_step(path_id(), b); }
        inline void set_prev_id(const uint64_t& i) { data[1] = i; }
        inline void set_prev_rank(const uint64_t& i) { data[2] = i; }
        inline void set_next_id(const uint64_t& i) { data[3] = i; }
        inline void set_next_rank(const uint64_t& i) { data[4] = i; }
    };
    const layout_t get_layout(void) const;
    const layout_t set_layout(layout_t layout);
    void init(const std::string& sequence);
    uint64_t sequence_size(void) const;
    const std::string sequence(void) const;
    void set_sequence(const std::string& seq);
    uint64_t edge_count_offset(void) const;
    uint64_t edge_count(void) const;
    uint64_t edge_offset(void) const;
    const std::vector<uint64_t> edges(void) const;
    void add_edge(const uint64_t& relative_id, const uint64_t& edge_type);
    void remove_edge(const uint64_t& rank);
    void set_edge_count(const uint64_t& count, const layout_t& layout);
    void add_path_step(const uint64_t& path_id, const bool& is_rev,
                       const uint64_t& prev_id, const uint64_t& prev_rank,
                       const uint64_t& next_id, const uint64_t& next_rank);
    void add_path_step(const step_t& step);
    void set_path_step(const uint64_t& rank, const uint64_t& path_id, const bool& is_rev,
                       const uint64_t& prev_id, const uint64_t& prev_rank,
                       const uint64_t& next_id, const uint64_t& next_rank);
    void set_path_step(const uint64_t& rank, const step_t& step);
    void flip_paths(void);
    const std::vector<node_t::step_t> get_path_steps(void) const;
    const step_t get_path_step(const uint64_t& rank) const;
    void remove_path_step(const uint64_t& rank);
    void clear_path_steps(void);
    uint64_t path_count(void) const;
    uint64_t serialize(std::ostream& out) const;
    void load(std::istream& in);
    void display(void) const;

    // path step helpers
    inline static uint64_t pack_step(const uint64_t& path_id, const bool& is_rev) {
        assert(path_id < (0x1ULL << 63));
        return ((path_id << 1) | (is_rev ? 1 : 0));
    }
    inline static uint64_t step_path_id(const uint64_t& step) {
        return step >> 1;
    }
    inline static bool step_is_rev(const uint64_t& step) {
        return step & 1;
    }

    // Some offset ints used in edge storage
    /*
    const static uint64_t EDGE_RECORD_LENGTH = 2;
    const static uint64_t EDGE_DELTA_OFFSET = 0;
    const static uint64_t EDGE_TYPE_OFFSET = 1;
    */
};

}

#endif
