#ifndef dgraph_node_hpp
#define dgraph_node_hpp

#include <iostream>
#include <cstdint>
#include <string>
#include <handlegraph/util.hpp>
#include "varint.hpp"

namespace dsgvg {

using namespace handlegraph;
using id_t = handlegraph::id_t;

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
    node_t(id_t id, const std::string& sequence) {
        // record the id and sequence size
        varint::encode({(uint64_t)id, (uint64_t)sequence.size()}, bytes);
        // and store the sequence
        bytes.insert(bytes.end(), sequence.begin(), sequence.end());
        varint::encode({0}, bytes); // record that there are no edges
        varint::encode({0}, bytes); // record that there are no paths
        //shrink_to_fit(); // TODO checkme
        //std::cerr << "created " << id << " " << sequence << std::endl;
        //display();
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
        uint64_t id = 0;
        uint64_t seq_start = 0;
        uint64_t seq_len = 0;
        uint64_t edge_start = 0;
        uint64_t edge_count = 0;
        uint64_t edge_bytes = 0;
        uint64_t path_start = 0;
        uint64_t path_count = 0;
        uint64_t path_bytes = 0;
    };
    struct step_t {
        uint64_t path_id = 0;
        bool is_rev = false;
        uint64_t prev_id = 0;
        uint64_t prev_rank = 0;
        uint64_t next_id = 0;
        uint64_t next_rank = 0;
    };
    const uint8_t EDGE_RECORD_LENGTH = 2;
    const uint8_t PATH_RECORD_LENGTH = 5;
    const layout_t get_layout(void) const;
    const layout_t get_seq_layout(void) const;
    const layout_t get_seq_edge_layout(void) const;
    const layout_t get_seq_edge_path_layout(void) const;
    id_t id(void) const;
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
    void set_path_step(const uint64_t& rank, const uint64_t& path_id, const bool& is_rev,
                       const uint64_t& prev_id, const uint64_t& prev_rank,
                       const uint64_t& next_id, const uint64_t& next_rank);
    void flip_paths(void);
    //step_t get_path_step_id(const uint64_t& rank) const;
    //step_t get_path_step_rev(const uint64_t& rank) const;
    const std::vector<node_t::step_t> get_path_steps(void) const;
    const step_t get_path_step(const uint64_t& rank) const;
    void remove_path_step(const uint64_t& rank);
    void clear_path_steps(void);
    void set_path_count(uint64_t count, const layout_t& layout);
    uint64_t path_count(void) const;
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
