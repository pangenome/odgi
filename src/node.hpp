#ifndef dgraph_node_hpp
#define dgraph_node_hpp

#include <iostream>
#include <cstdint>
#include <string>
#include "handle.hpp"
#include "varint.hpp"

namespace dsgvg {

/// A node object with the sequence, its edge lists, and paths
struct node_t : std::vector<uint8_t> {
    node_t(void);
    node_t(id_t id, const std::string& sequence);
    struct layout_t {
        uint64_t id = 0;
        uint64_t seq_start = 0;
        uint64_t seq_len = 0;
        uint64_t edge_start = 0;
        uint64_t edge_count = 0;
        // uint64_t path_start = 0; // TODO
    };
    const layout_t get_layout(void) const;
    const layout_t get_seq_layout(void) const;
    // stored uncompressed for performance
    id_t id(void) const;
    uint64_t sequence_size(void) const;
    const std::string sequence(void) const;
    void set_sequence(const std::string& seq);
    uint64_t edge_count_offset(void) const;
    uint64_t edge_count(void) const;
    uint64_t edge_offset(void) const;
    const std::vector<uint64_t> edges(void) const;
    void add_edge(uint64_t relative_id, uint64_t edge_type);
    void remove_edge(uint64_t i);
    void set_edge_count(uint64_t count, const layout_t& layout);
    void display(void) const;

    // Some offset ints used in 
    const static uint64_t EDGE_RECORD_LENGTH = 2;
    const static uint64_t EDGE_DELTA_OFFSET = 0;
    const static uint64_t EDGE_TYPE_OFFSET = 1;
};

}

#endif
