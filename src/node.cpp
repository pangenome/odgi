#include "node.hpp"
#include "varint.hpp"
#include <cassert>

namespace odgi {

uint64_t node_t::sequence_size(void) const {
    return sequence.size();
}

const std::string node_t::get_sequence(void) const {
    return sequence;
}

void node_t::set_sequence(const std::string& seq) {
    sequence = seq;
}

const dyn::hacked_vector& node_t::get_edges(void) const {
    return edges;
}

void node_t::add_edge(const uint64_t& relative_id, const uint64_t& edge_type) {
    //std::cerr << "add edge " << "relative_id " << relative_id << " edge_type " << edge_type << std::endl;
    edges.push_back(relative_id);
    edges.push_back(edge_type);
}

void node_t::remove_edge(const uint64_t& rank) {
    assert(rank < edge_count());
    uint64_t offset = EDGE_RECORD_LENGTH*rank;
    for (uint8_t i = 0; i < EDGE_RECORD_LENGTH; ++i) {
        edges.remove(offset);
    }
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
    sequence.clear();
    clear_edges();
    clear_path_steps();
}

void node_t::clear_edges(void) {
    dyn::hacked_vector null_iv;
    edges = null_iv;
}

void node_t::clear_path_steps(void) {
    dyn::hacked_vector null_iv;
    path_steps = null_iv;
}

uint64_t node_t::serialize(std::ostream& out) const {
    uint64_t written = 0;
    size_t seq_size = sequence.size();
    out.write((char*)&seq_size, sizeof(size_t));
    written += sizeof(size_t);
    out << sequence;
    written += sequence.size();
    written += edges.serialize(out);
    written += path_steps.serialize(out);
    return written;
}

void node_t::load(std::istream& in) {
    size_t seq_size;
    in.read((char*)&seq_size, sizeof(size_t));
    sequence.resize(seq_size);
    in.read((char*)sequence.c_str(), seq_size);
    edges.load(in);
    path_steps.load(in);
}

void node_t::display(void) const {
    std::cerr << "seq " << sequence << " "
              << "edge_count " << edge_count() << " "
              << "path_count " << path_count() << " | ";
    if (edge_count()) {
        for (uint64_t i = 0; i < edge_count(); ++i) {
            std::cerr
                << edges.at(i) << ":"
                << edges.at(i+1) << " ";
        }
    }
    std::cerr << "| ";
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
