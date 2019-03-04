#include "node.hpp"

namespace odgi {

const node_t::layout_t node_t::get_layout(void) const {
    node_t::layout_t layout;
    varint::decode((uint8_t*)bytes.data(), &layout.data[0], 5);
    return layout;
}

const node_t::layout_t node_t::set_layout(node_t::layout_t layout) {
    // how big is our layout?
    uint64_t old_size = varint::bytes(bytes.data(), 5);
    uint64_t new_size = varint::length(&layout.data[0], 5);
    if (new_size < old_size) {
        bytes.erase(bytes.begin(), bytes.begin()+(old_size-new_size));
    } else if (new_size > old_size) {
        bytes.insert(bytes.begin(), new_size-old_size, 0);
    }
    layout.set_layout_bytes(new_size);
    varint::encode(layout.data, bytes.data(), 5);
    return layout;
}

void node_t::init(const std::string& sequence) {
    set_sequence(sequence);
}

uint64_t node_t::sequence_size(void) const {
    node_t::layout_t layout = get_layout();
    return layout.seq_bytes();
}

const std::string node_t::sequence(void) const {
    node_t::layout_t layout = get_layout();
    const std::string res((char*)bytes.data()+layout.layout_bytes()+layout.seq_bytes(), layout.seq_bytes());
    return res;
}

void node_t::set_sequence(const std::string& seq) {
    node_t::layout_t layout = get_layout();
    if (seq.size() > layout.seq_bytes()) {
        bytes.insert(bytes.begin()+layout.seq_start(), seq.size() - layout.seq_bytes(), 0);
        layout.set_seq_bytes(seq.size());
        layout = set_layout(layout);
    } else if (seq.size() < layout.seq_bytes()) {
        bytes.erase(bytes.begin()+layout.seq_start(), bytes.begin()+layout.seq_start()+(layout.seq_bytes()-seq.size()));;
        layout.set_seq_bytes(seq.size());
        layout = set_layout(layout);
    }
    memcpy(bytes.data()+layout.seq_start(), seq.c_str(), seq.size());
}

const std::vector<uint64_t> node_t::edges(void) const {
    std::vector<uint64_t> res;
    const node_t::layout_t& layout = get_layout();
    if (layout.edge_count()) {
        res.resize(layout.edge_count()*EDGE_RECORD_LENGTH);
        varint::decode((uint8_t*)bytes.data()+layout.edge_start(), res.data(), layout.edge_count()*EDGE_RECORD_LENGTH);
    }
    assert(res.size() == layout.edge_count);
    return res;
}

void node_t::add_edge(const uint64_t& relative_id, const uint64_t& edge_type) {
    node_t::layout_t layout = get_layout();
    uint64_t edge_bytes = varint::length({relative_id, edge_type});
    bytes.insert(bytes.begin()+layout.edge_start(), edge_bytes, 0);
    varint::encode({relative_id, edge_type}, bytes.data()+layout.edge_start());
    layout.set_edge_bytes(layout.edge_bytes() + edge_bytes);
    layout.set_edge_count(layout.edge_count() + 1);
    set_layout(layout);
}

void node_t::remove_edge(const uint64_t& rank) {
    assert(rank < edge_count());
    node_t::layout_t layout = get_layout();
    if (rank > layout.edge_count()) assert(false);
    uint64_t edge_offset = layout.edge_start() + varint::bytes(bytes.data()+layout.edge_start(), EDGE_RECORD_LENGTH*(rank-1));
    // a bit redundant
    uint64_t j = varint::bytes(bytes.data()+edge_offset, EDGE_RECORD_LENGTH);
    bytes.erase(bytes.begin()+edge_offset, bytes.begin()+edge_offset+j);
    layout.set_edge_count(layout.edge_count()-1);
    layout.set_edge_bytes(layout.edge_bytes()-j);
    set_layout(layout);
}

uint64_t node_t::edge_count(void) const {
    return get_layout().edge_count();
}

void node_t::add_path_step(const uint64_t& path_id, const bool& is_rev,
                           const uint64_t& prev_id, const uint64_t& prev_rank,
                           const uint64_t& next_id, const uint64_t& next_rank) {
    node_t::step_t step;
    step.data[0] = pack_step(path_id, is_rev);
    step.data[1] = prev_id;
    step.data[2] = prev_rank;
    step.data[3] = next_id;
    step.data[4] = next_rank;
    add_path_step(step);
}

void node_t::add_path_step(const node_t::step_t& step) {
    layout_t layout = get_layout();
    layout.set_path_count(layout.path_count()+1);
    layout = set_layout(layout);
    varint::encode(step.data, bytes.data(), 5);
}

const std::vector<node_t::step_t> node_t::get_path_steps(void) const {
    layout_t layout = get_layout();
    if (layout.path_count() == 0) return {};
    std::vector<node_t::step_t> steps(layout.path_count());
    uint8_t* target = (uint8_t*)bytes.data()+layout.path_start();
    for (uint64_t i = 0; i < layout.path_count(); ++i) {
        auto& step = steps[i];
        target = varint::decode(target, &step.data[0], PATH_RECORD_LENGTH);
    }
    return steps;
}

const node_t::step_t node_t::get_path_step(const uint64_t& rank) const {
    layout_t layout = get_layout();
    node_t::step_t step;
    varint::decode(varint::seek((uint8_t*)bytes.data() + layout.path_start(), PATH_RECORD_LENGTH*rank),
                   &step.data[0], PATH_RECORD_LENGTH);
    return step;
}

void node_t::set_path_step(const uint64_t& rank, const step_t& step) {
    layout_t layout = get_layout();
    uint8_t* target = varint::seek(bytes.data() + layout.path_start(), PATH_RECORD_LENGTH*rank);
    uint64_t old_size = varint::bytes(target, PATH_RECORD_LENGTH);
    uint64_t new_size = varint::length(&step.data[0], PATH_RECORD_LENGTH);
    if (new_size > old_size) {
        // insert
        std::vector<uint8_t>::iterator it(target);
        bytes.insert(it, new_size - old_size, 0);
    } else if (new_size < old_size) {
        std::vector<uint8_t>::iterator it(target);
        bytes.erase(it, it + (old_size - new_size));
    }
    varint::encode(step.data, target, 5);
}

void node_t::flip_paths(void) {
    const std::vector<node_t::step_t> steps = get_path_steps();
    // remove all path steps
    node_t::layout_t layout = get_layout();
    bytes.erase(bytes.begin()+layout.path_start(), bytes.end());
    // flip them and replace
    for (auto& step : steps) {
        bytes.resize(bytes.size()+5);
        varint::encode(step.data, (uint8_t*)bytes.data()+bytes.size(), 5);
    }
}

void node_t::remove_path_step(const uint64_t& rank) {
    node_t::layout_t layout = get_layout();
    if (rank > layout.path_count()) assert(false);
    uint8_t* i = varint::seek(bytes.data()+layout.path_start(), PATH_RECORD_LENGTH*(rank-1));
    uint8_t* j = varint::seek(i, PATH_RECORD_LENGTH);
    bytes.erase((std::vector<uint8_t>::iterator)i,
                (std::vector<uint8_t>::iterator)j);
    layout.set_path_count(layout.path_count()-1);
    set_layout(layout);
}

void node_t::clear_path_steps(void) {
    node_t::layout_t layout = get_layout();
    bytes.erase(bytes.begin()+layout.path_start(), bytes.end());
    layout.set_path_count(0);
    set_layout(layout);
}

uint64_t node_t::path_count(void) const {
    layout_t layout = get_layout();
    return layout.path_count();
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
    layout_t layout = get_layout();
    std::cerr << layout.layout_bytes() << " "
              << layout.seq_bytes() << " "
              << layout.edge_start() << " "
              << layout.edge_count() << " "
              << layout.edge_bytes() << " "
              << layout.path_start() << " "
              << layout.path_count() << " | ";
    for (auto i : bytes) {
        std::cerr << (int) i << " ";
    }
    std::cerr << std::endl;
}

}
