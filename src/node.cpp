#include "node.hpp"

namespace odgi {

const node_t::layout_t node_t::get_seq_layout(void) const {
    node_t::layout_t layout;
    layout.id = id();
    layout.seq_start = varint::length(layout.id);
    varint::decode(bytes.data()+layout.seq_start, &layout.seq_len, 1);
    return layout;
}

const node_t::layout_t node_t::get_seq_edge_layout(void) const {
    node_t::layout_t layout;
    layout.id = id();
    layout.seq_start = varint::length(layout.id);
    varint::decode(bytes.data()+layout.seq_start, &layout.seq_len, 1);
    layout.edge_start = layout.seq_start + varint::length(layout.seq_len) + layout.seq_len;
    varint::decode(bytes.data()+layout.edge_start, &layout.edge_count, 1);
    layout.edge_bytes = varint::bytes(bytes.data()+layout.edge_start, 1+layout.edge_count*EDGE_RECORD_LENGTH);
    return layout;
}

const node_t::layout_t node_t::get_seq_edge_path_layout(void) const {
    node_t::layout_t layout;
    layout.id = id();
    layout.seq_start = varint::length(layout.id);
    varint::decode(bytes.data()+layout.seq_start, &layout.seq_len, 1);
    layout.edge_start = layout.seq_start + varint::length(layout.seq_len) + layout.seq_len;
    varint::decode(bytes.data()+layout.edge_start, &layout.edge_count, 1);
    layout.edge_bytes = varint::bytes(bytes.data()+layout.edge_start, 1+layout.edge_count*EDGE_RECORD_LENGTH);
    layout.path_start = layout.edge_start + layout.edge_bytes;
    varint::decode(bytes.data()+layout.path_start, &layout.path_count, 1);
    layout.path_bytes = varint::length(layout.path_count)+varint::bytes(bytes.data()+layout.path_start+varint::length(layout.path_count), layout.path_count*PATH_RECORD_LENGTH);
    return layout;
}

void node_t::init(const id_t& id, const std::string& sequence) {
    // record the id and sequence size
    varint::encode({(uint64_t)id, (uint64_t)sequence.size()}, bytes);
    // and store the sequence
    bytes.insert(bytes.end(), sequence.begin(), sequence.end());
    varint::encode({0}, bytes); // record that there are no edges
    varint::encode({0}, bytes); // record that there are no paths
}

id_t node_t::id(void) const {
    uint64_t res;
    varint::decode(bytes.data(), &res, 1);
    return (id_t)res;
}

uint64_t node_t::sequence_size(void) const {
    uint64_t len = 0;
    varint::decode(bytes.data()+varint::length(id()), &len, 1);
    return len;
}

const std::string node_t::sequence(void) const {
    node_t::layout_t layout = get_seq_layout();
    const std::string res((char*)bytes.data()+layout.seq_start+varint::length(layout.seq_len), layout.seq_len);
    return res;
}

void node_t::set_sequence(const std::string& seq) {
    node_t::layout_t layout = get_seq_layout();
    bytes.erase(bytes.begin()+layout.seq_start, bytes.begin()+layout.seq_start+varint::length(layout.seq_len)+layout.seq_len);
    std::vector<uint8_t> seq_record = varint::encode({seq.size()});
    seq_record.reserve(seq_record.size()+seq.size());
    for (auto c : seq) seq_record.push_back(c);
    bytes.insert(bytes.begin()+layout.seq_start, seq_record.begin(), seq_record.end());
}

const std::vector<uint64_t> node_t::edges(void) const {
    std::vector<uint64_t> res;
    node_t::layout_t layout = get_seq_edge_layout();
    if (layout.edge_count) {
        res.resize(layout.edge_count*EDGE_RECORD_LENGTH);
        varint::decode(bytes.data()+layout.edge_start+varint::length(layout.edge_count), res.data(), layout.edge_count*EDGE_RECORD_LENGTH);
    }
    assert(res.size() == layout.edge_count);
    return res;
}

void node_t::add_edge(const uint64_t& relative_id, const uint64_t& edge_type) {
    node_t::layout_t layout = get_seq_edge_path_layout();
    std::vector<uint8_t> edge = varint::encode({relative_id, edge_type});
    bytes.insert(bytes.begin()+layout.edge_start+varint::length(layout.edge_count), edge.begin(), edge.end());
    set_edge_count(layout.edge_count+1, layout);
}

void node_t::remove_edge(const uint64_t& rank) {
    assert(rank < edge_count());
    node_t::layout_t layout = get_seq_edge_layout();
    if (rank > layout.edge_count) assert(false);
    uint64_t edge_offset = layout.edge_start + varint::length(layout.edge_count);
    
    uint64_t i = varint::bytes(bytes.data()+edge_offset, EDGE_RECORD_LENGTH*rank);
    // a bit redundant
    uint64_t j = i+varint::bytes(bytes.data()+edge_offset+i, EDGE_RECORD_LENGTH);
    bytes.erase(bytes.begin()+edge_offset+i, bytes.begin()+edge_offset+j);
    set_edge_count(layout.edge_count-1, layout);
}

void node_t::set_edge_count(const uint64_t& count, const node_t::layout_t& layout) {
    uint8_t count_encoding_delta = varint::length(count) - varint::length(layout.edge_count);
    // check if we need to realloc the edge count storage
    if (count_encoding_delta < 0) {
        bytes.erase(bytes.begin()+layout.edge_start, bytes.begin()+layout.edge_start+varint::length(layout.edge_count));
        std::vector<uint8_t> x(count);
        bytes.insert(bytes.begin()+layout.edge_start, x.begin(), x.end());
    } else if (count_encoding_delta > 0) {
        std::vector<uint8_t> x(count_encoding_delta);
        bytes.insert(bytes.begin()+layout.edge_start, x.begin(), x.end());
    }
    std::vector<uint8_t> c = varint::encode({count});
    memcpy(bytes.data()+layout.edge_start, c.data(), c.size());
}

uint64_t node_t::edge_count(void) const {
    layout_t layout = get_seq_edge_layout();
    return layout.edge_count;
}

void node_t::add_path_step(const uint64_t& path_id, const bool& is_rev,
                           const uint64_t& prev_id, const uint64_t& prev_rank,
                           const uint64_t& next_id, const uint64_t& next_rank) {
    layout_t layout = get_seq_edge_path_layout();
    std::vector<uint8_t> step = varint::encode({pack_step(path_id, is_rev), prev_id, prev_rank, next_id, next_rank});
    bytes.insert(bytes.begin()+layout.path_start+layout.path_bytes, step.begin(), step.end());
    set_path_count(layout.path_count+1, layout);
    //display();
}

const std::vector<node_t::step_t> node_t::get_path_steps(void) const {
    layout_t layout = get_seq_edge_path_layout();
    std::vector<uint64_t> step_v;
    if (layout.path_count == 0) return {};
    step_v.resize(PATH_RECORD_LENGTH*layout.path_count);
    varint::decode(bytes.data()+layout.path_start+varint::length(layout.path_count), step_v.data(), PATH_RECORD_LENGTH*layout.path_count);
    assert(step_v.size() == layout.path_count*PATH_RECORD_LENGTH);
    std::vector<node_t::step_t> steps;
    for (uint64_t i = 0; i < layout.path_count*PATH_RECORD_LENGTH; i+=PATH_RECORD_LENGTH) {
        //auto& step = steps[j++];
        steps.emplace_back();
        auto& step = steps.back();
        step.path_id = step_path_id(step_v.at(i));
        step.is_rev = step_is_rev(step_v.at(i));
        step.prev_id = step_v.at(i+1);
        step.prev_rank = step_v.at(i+2);
        step.next_id = step_v.at(i+3);
        step.next_rank = step_v.at(i+4);
    }
    return steps;
}

const node_t::step_t node_t::get_path_step(const uint64_t& rank) const {
    layout_t layout = get_seq_edge_path_layout();
    std::vector<uint64_t> step_v;
    step_v.resize(PATH_RECORD_LENGTH);
    varint::decode(bytes.data()+layout.path_start+varint::bytes(bytes.data()+layout.path_start, 1+PATH_RECORD_LENGTH*rank), step_v.data(), PATH_RECORD_LENGTH);
    node_t::step_t step;
    step.path_id = step_path_id(step_v.at(0));
    step.is_rev = step_is_rev(step_v.at(0));
    step.prev_id = step_v.at(1);
    step.prev_rank = step_v.at(2);
    step.next_id = step_v.at(3);
    step.next_rank = step_v.at(4);
    return step;
}

void node_t::set_path_step(const uint64_t& rank, const uint64_t& path_id, const bool& is_rev,
                           const uint64_t& prev_id, const uint64_t& prev_rank,
                           const uint64_t& next_id, const uint64_t& next_rank) {
    layout_t layout = get_seq_edge_path_layout();
    // erase the step
    std::vector<uint64_t> step_v;
    uint64_t step_begin = layout.path_start+varint::bytes(bytes.data()+layout.path_start, 1+PATH_RECORD_LENGTH*(rank));
    uint64_t step_end = layout.path_start+varint::bytes(bytes.data()+layout.path_start, 1+PATH_RECORD_LENGTH*(rank+1));
    bytes.erase(bytes.begin()+step_begin, bytes.begin()+step_end);
    // insert the step
    std::vector<uint8_t> step = varint::encode({pack_step(path_id, is_rev), prev_id, prev_rank, next_id, next_rank});
    bytes.insert(bytes.begin()+step_begin, step.begin(), step.end());
}

void node_t::flip_paths(void) {
    const std::vector<node_t::step_t> steps = get_path_steps();
    // remove all path steps
    node_t::layout_t layout = get_seq_edge_path_layout();
    bytes.erase(bytes.begin()+layout.path_start+varint::length(layout.path_count), bytes.end());
    // flip them and replace
    for (auto& step : steps) {
        //std::vector<uint8_t> step_v = varint::encode({pack_step(step.path_id, !step.is_rev), step.prev_id, step.prev_rank, step.next_id, step.next_rank});
        //bytes.insert(bytes.end(), step_v.bytes.begin(), step_v.bytes.end());
        varint::encode({pack_step(step.path_id, !step.is_rev), step.prev_id, step.prev_rank, step.next_id, step.next_rank}, bytes);
    }
}

void node_t::remove_path_step(const uint64_t& rank) {
    node_t::layout_t layout = get_seq_edge_path_layout();
    if (rank > layout.path_count) assert(false);
    uint64_t path_offset = layout.path_start + varint::length(layout.path_count);
    uint64_t i = varint::bytes(bytes.data()+path_offset, PATH_RECORD_LENGTH*(rank));
    // a bit redundant
    uint64_t j = varint::bytes(bytes.data()+path_offset, PATH_RECORD_LENGTH*(rank+1));
    bytes.erase(bytes.begin()+path_offset+i, bytes.begin()+path_offset+j);
    set_path_count(layout.path_count-1, layout);
}

void node_t::clear_path_steps(void) {
    node_t::layout_t layout = get_seq_edge_path_layout();
    bytes.erase(bytes.begin()+layout.path_start, bytes.end());
    // set counter for 0 paths
    varint::encode({0}, bytes);
}

void node_t::set_path_count(uint64_t count, const layout_t& layout) {
    uint8_t count_encoding_delta = varint::length(count) - varint::length(layout.path_count);
    // check if we need to realloc the edge count storage
    if (count_encoding_delta < 0) {
        bytes.erase(bytes.begin()+layout.path_start, bytes.begin()+layout.path_start+varint::length(layout.path_count));
        std::vector<uint8_t> x(count);
        bytes.insert(bytes.begin()+layout.path_start, x.begin(), x.end());
    } else if (count_encoding_delta > 0) {
        std::vector<uint8_t> x(count_encoding_delta);
        bytes.insert(bytes.begin()+layout.path_start, x.begin(), x.end());
    }
    std::vector<uint8_t> c = varint::encode({count});
    memcpy(bytes.data()+layout.path_start, c.data(), c.size());
}

uint64_t node_t::path_count(void) const {
    layout_t layout = get_seq_edge_path_layout();
    return layout.path_count;
}

uint64_t node_t::serialize(std::ostream& out) const {
    uint64_t written = 0;
    uint64_t node_size = bytes.size();
    out.write((char*)&node_size, sizeof(node_size));
    written += sizeof(uint64_t);
    out.write((char*)bytes.data(), node_size*sizeof(uint8_t));
    written += sizeof(uint8_t)*node_size;
    return written;
}

void node_t::load(std::istream& in) {
    uint64_t node_size = 0;
    in.read((char*)&node_size, sizeof(node_size));
    bytes.resize(node_size);
    in.read((char*)bytes.data(), node_size*sizeof(uint8_t));
}

void node_t::display(void) const {
    layout_t layout = get_seq_edge_path_layout();
    std::cerr << layout.id << " "
              << layout.seq_start << " "
              << layout.seq_len << " "
              << layout.edge_start << " "
              << layout.edge_count << " "
              << layout.edge_bytes << " "
              << layout.path_start << " "
              << layout.path_count << " "
              << layout.path_bytes << " | ";
    for (auto i : bytes) {
        std::cerr << (int) i << " ";
    }
    std::cerr << std::endl;
}

}
