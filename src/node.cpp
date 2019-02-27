#include "node.hpp"

namespace dsgvg {

node_t::node_t(void) { }

node_t::node_t(id_t id, const std::string& sequence) {
    // record the id and sequence size
    varint::encode({(uint64_t)id, (uint64_t)sequence.size()}, *this);
    // and store the sequence
    insert(end(), sequence.begin(), sequence.end());
    varint::encode({0}, *this); // record that there are no edges
}

const node_t::layout_t node_t::get_layout(void) const {
    node_t::layout_t layout;
    layout.id = id();
    layout.seq_start = varint::length(layout.id);
    varint::decode(data()+layout.seq_start, &layout.seq_len, 1);
    layout.edge_start = layout.seq_start + varint::length(layout.seq_len) + layout.seq_len;
    varint::decode(data()+layout.edge_start, &layout.edge_count, 1);
    return layout;
}

const node_t::layout_t node_t::get_seq_layout(void) const {
    node_t::layout_t layout;
    layout.id = id();
    layout.seq_start = varint::length(layout.id);
    varint::decode(data()+layout.seq_start, &layout.seq_len, 1);
    return layout;
}

id_t node_t::id(void) const {
    uint64_t res;
    varint::decode(data(), &res, 1);
    return (id_t)res;
}

uint64_t node_t::sequence_size(void) const {
    uint64_t len = 0;
    varint::decode(data()+varint::length(id()), &len, 1);
    return len;
}

const std::string node_t::sequence(void) const {
    node_t::layout_t layout = get_seq_layout();
    const std::string res((char*)data()+layout.seq_start+varint::length(layout.seq_len), layout.seq_len);
    return res;
}

void node_t::set_sequence(const std::string& seq) {
    node_t::layout_t layout = get_seq_layout();
    erase(begin()+layout.seq_start, begin()+layout.seq_start+varint::length(layout.seq_len)+layout.seq_len);
    std::vector<uint8_t> seq_record = varint::encode({seq.size()});
    seq_record.reserve(seq_record.size()+seq.size());
    for (auto c : seq) seq_record.push_back(c);
    insert(begin()+layout.seq_start, seq_record.begin(), seq_record.end());
}

const std::vector<uint64_t> node_t::edges(void) const {
    std::vector<uint64_t> res;
    node_t::layout_t layout = get_layout();
    if (layout.edge_count) {
        res.resize(layout.edge_count*2);
        varint::decode(data()+layout.edge_start+varint::length(layout.edge_count), res.data(), layout.edge_count*2);
    }
    return res;
}

void node_t::add_edge(uint64_t relative_id, uint64_t edge_type) {
    node_t::layout_t layout = get_layout();
    set_edge_count(layout.edge_count+1, layout);
    varint::encode({relative_id, edge_type}, *this);
}

void node_t::remove_edge(uint64_t i) {
    node_t::layout_t layout = get_layout();
    if (i > layout.edge_count*2) assert(false);
    uint64_t edge_offset = layout.edge_start + varint::length(layout.edge_count);
    erase(begin()+edge_offset+i, begin()+edge_offset+i+2);
    set_edge_count(layout.edge_count-1, layout);
}

void node_t::set_edge_count(uint64_t count, const node_t::layout_t& layout) {
    uint8_t count_encoding_delta = varint::length(count) - varint::length(layout.edge_count);
    // check if we need to realloc the edge count storage
    if (count_encoding_delta < 0) {
        erase(begin()+layout.edge_start, begin()+layout.edge_start+varint::length(layout.edge_count));
        std::vector<uint8_t> x(count);
        insert(begin()+layout.edge_start, x.begin(), x.end());
    } else if (count_encoding_delta > 0) {
        std::vector<uint8_t> x(count_encoding_delta);
        insert(begin()+layout.edge_start, x.begin(), x.end());
    }
    std::vector<uint8_t> c = varint::encode({count});
    memcpy(data()+layout.edge_start, c.data(), c.size());
}

void node_t::display(void) const {
    layout_t layout = get_layout();
    std::cerr << layout.id << " "
              << layout.seq_start << " "
              << layout.seq_len << " "
              << layout.edge_start << " "
              << layout.edge_count << " |";
    for (auto i : *this) {
        std::cerr << (int) i << " ";
    }
    std::cerr << std::endl;
}

}
