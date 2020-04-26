#include "node.hpp"
#include "varint.hpp"
#include <cassert>

namespace odgi {

uint64_t node_t::sequence_size(void) const {
    return seq_bytes();
}

const std::string node_t::sequence(void) const {
    const std::string res((char*)bytes.data()+seq_start(), seq_bytes());
    return res;
}

void node_t::set_sequence(const std::string& seq) {
    if (seq.size() > seq_bytes()) {
        bytes.reserve(bytes.size()+seq.size()-seq_bytes());
        bytes.insert(bytes.begin()+seq_start(), seq.size() - seq_bytes(), 0);
        set_seq_bytes(seq.size());
    } else if (seq.size() < seq_bytes()) {
        bytes.erase(bytes.begin()+seq_start(), bytes.begin()+seq_start()+(seq_bytes()-seq.size()));;
        set_seq_bytes(seq.size());
    }
    memcpy(bytes.data()+seq_start(), seq.c_str(), seq.size());
}

std::vector<uint64_t> node_t::edges(void) const {
    std::vector<uint64_t> res;
    if (edge_count()) {
        res.resize(edge_count()*EDGE_RECORD_LENGTH);
        sqvarint::decode(res.data(),
                       (uint8_t*)bytes.data()+edge_start(),
                       edge_count()*EDGE_RECORD_LENGTH);
    }
    return res;
}

void node_t::add_edge(const uint64_t& relative_id, const uint64_t& edge_type) {
    //std::cerr << "add edge " << "relative_id " << relative_id << " edge_type " << edge_type << std::endl;
    uint64_t add_edge_bytes = sqvarint::length({relative_id, edge_type});
    bytes.reserve(bytes.size()+add_edge_bytes);
    bytes.insert(bytes.begin()+edge_start(), add_edge_bytes, 0);
    sqvarint::encode({relative_id, edge_type}, bytes.data()+edge_start());
    set_edge_bytes(edge_bytes() + add_edge_bytes);
    set_edge_count(edge_count() + 1);
}

void node_t::remove_edge(const uint64_t& rank) {
    assert(rank < edge_count());
    uint64_t edge_offset = edge_start() + sqvarint::bytes(bytes.data()+edge_start(), EDGE_RECORD_LENGTH*rank);
    // a bit redundant
    uint64_t j = sqvarint::bytes(bytes.data()+edge_offset, EDGE_RECORD_LENGTH);
    bytes.erase(bytes.begin()+edge_offset, bytes.begin()+edge_offset+j);
    set_edge_count(edge_count()-1);
    set_edge_bytes(edge_bytes()-j);
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
    for (uint8_t i = 0; i < PATH_RECORD_LENGTH; ++i) {
        path_steps.push_back(step.data[i]);
    }
}

const std::vector<node_t::step_t> node_t::get_path_steps(void) const {
    uint64_t n_paths = path_count();
    if (n_paths == 0) return {};
    std::vector<node_t::step_t> steps(n_paths);
    for (uint64_t i = 0; i < n_paths; ++i) {
        steps[i] = get_path_step(i);
    }
    return steps;
}

const node_t::step_t node_t::get_path_step(const uint64_t& rank) const {
    if (rank >= path_count()) assert(false);
    node_t::step_t step;
    uint64_t offset = PATH_RECORD_LENGTH*rank;
    for (uint8_t i = 0; i < PATH_RECORD_LENGTH; ++i) {
        step.data[i] = path_steps.at(offset+i);
    }
    return step;
}

void node_t::set_path_step(const uint64_t& rank, const step_t& step) {
    if (rank >= path_count()) assert(false);
    uint64_t offset = PATH_RECORD_LENGTH*rank;
    for (uint8_t i = 0; i < PATH_RECORD_LENGTH; ++i) {
        path_steps[offset+i] = step.data[i];
    }
}

std::pair<std::map<uint64_t, std::pair<uint64_t, bool>>, // path fronts
          std::map<uint64_t, std::pair<uint64_t, bool>>> // path backs
node_t::flip_paths(const uint64_t& start_marker,
                   const uint64_t& end_marker) {
    const std::vector<node_t::step_t> steps = get_path_steps();
    // remove all path steps
    clear_path_steps();
    // flip them and replace, recording which path starts and ends should be rewritten
    std::pair<std::map<uint64_t, std::pair<uint64_t, bool>>,
              std::map<uint64_t, std::pair<uint64_t, bool>>> path_start_end_rewrites;
    uint64_t rank = 0;
    for (auto& step : steps) {
        // flip the step
        step_t flipped(step.path_id(), !step.is_rev(),
                       step.prev_id(), step.prev_rank(),
                       step.next_id(), step.next_rank());
        if (step.prev_id() == start_marker) {
            path_start_end_rewrites.first[step.path_id()] = std::make_pair(rank, !step.is_rev());
        }
        if (step.next_id() == end_marker) {
            path_start_end_rewrites.second[step.path_id()] = std::make_pair(rank, !step.is_rev());
        }
        add_path_step(flipped);
        ++rank;
    }
    return path_start_end_rewrites;
}

void node_t::remove_path_step(const uint64_t& rank) {
    if (rank >= path_count()) assert(false);
    uint64_t offset = PATH_RECORD_LENGTH*rank;
    for (uint8_t i = 0; i < PATH_RECORD_LENGTH; ++i) {
        path_steps.remove(offset);
    }
}

void node_t::clear(void) {
    set_seq_bytes(0);
    set_edge_bytes(0);
    set_edge_count(0);
    bytes.clear();
    clear_path_steps();
}

void node_t::clear_path_steps(void) {
    dyn::hacked_vector null_iv;
    path_steps = null_iv;
}

uint64_t node_t::serialize(std::ostream& out) const {
    uint64_t written = 0;
    out.write((char*)&_seq_bytes, sizeof(uint32_t));
    out.write((char*)&_edge_bytes, sizeof(uint32_t));
    out.write((char*)&_edge_count, sizeof(uint32_t));
    written += sizeof(uint32_t)*4 + sizeof(uint8_t);
    uint64_t node_size = bytes.size();
    out.write((char*)&node_size, sizeof(node_size));
    written += sizeof(uint64_t);
    out.write((char*)bytes.data(), node_size*sizeof(uint8_t));
    written += sizeof(uint8_t)*node_size;
    written += path_steps.serialize(out);
    return written;
}

void node_t::load(std::istream& in) {
    in.read((char*)&_seq_bytes, sizeof(uint32_t));
    in.read((char*)&_edge_bytes, sizeof(uint32_t));
    in.read((char*)&_edge_count, sizeof(uint32_t));
    uint64_t node_size = 0;
    in.read((char*)&node_size, sizeof(node_size));
    bytes.resize(node_size);
    in.read((char*)bytes.data(), node_size*sizeof(uint8_t));
    path_steps.load(in);
}

void node_t::display(void) const {
    std::cerr << "self_bytes " << bytes.size() << " "
              << "seq_bytes " << seq_bytes() << " "
              << "seq " << sequence() << " "
              << "edge_start " << edge_start() << " "
              << "edge_count " << edge_count() << " "
              << "edge_bytes " << edge_bytes() << " "
              << "path_count " << path_count() << " | ";
    for (auto i : bytes) {
        std::cerr << (int) i << " ";
    }
    std::cerr << " | ";
    if (path_count()) {
        for (uint64_t i = 0; i < path_count(); ++i) {
            std::cerr
                << path_steps.at(i) << ":"
                << path_steps.at(i+1) << ":"
                << path_steps.at(i+2) << ":"
                << path_steps.at(i+3) << ":"
                << path_steps.at(i+4) << " ";
        }
    }
    std::cerr << std::endl;
}

}
