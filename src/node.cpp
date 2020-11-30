#include "node.hpp"
#include "varint.hpp"
#include <cassert>

namespace odgi {

uint64_t node_t::sequence_size(void) const {
    return sequence.size();
}

void node_t::set_sequence(const std::string& seq) {
    sequence = seq;
}

const std::string& node_t::get_sequence(void) const {
    return sequence;
}

/*
uint64_t node_t::delta_to_id(const uint64_t& base, const uint64_t& delta) {
    if (delta == 1) {
        return base;
    } else if (delta % 2 == 0) {
        return base + delta/2;
    } else { //if (delta-1 % 2 == 0) {
        return base - (delta-1)/2;
    }
}
*/

/*
dyn::hacked_vector& node_t::get_edges(void) const {
    return edges;
}
*/

void node_t::for_each_edge(const std::function<bool(nid_t other_id,
                                                    bool other_rev,
                                                    bool to_curr,
                                                    bool on_rev)>& func) const {
    for (uint64_t i = 0; i < edges.size(); ) {
        uint64_t other_id = edges.at(i++);
        uint8_t packed_edge = edges.at(i++);
        if (!func(other_id,
                  edge_helper::unpack_other_rev(packed_edge),
                  edge_helper::unpack_to_curr(packed_edge),
                  edge_helper::unpack_on_rev(packed_edge))) {
            break;
        }
    }
}

void node_t::add_edge(const uint64_t& other_id,
                      const bool& other_rev,
                      const bool& to_curr,
                      const bool& on_rev) {
    //std::cerr << "add edge " << "relative_id " << relative_id << " edge_type " << edge_type << std::endl;
    auto edge_type = edge_helper::pack(other_rev, to_curr, on_rev);
    edges.push_back(other_id);
    edges.push_back(edge_type);
}

bool node_t::remove_edge(const uint64_t& target_id,
                         const bool& target_rev,
                         const bool& ends_here,
                         const bool& is_rev) {
    std::cerr << "removing edge " << target_id << ":" << target_rev << ":" << ends_here << ":" << is_rev << std::endl;
    display();
    for (uint64_t i = 0; i < edges.size(); i+=EDGE_RECORD_LENGTH) {
        uint64_t other_id = edges.at(i);
        std::cerr << "other_id = " << other_id << std::endl;
        if (other_id == target_id) {
            uint8_t packed_edge = edges.at(i+1);
            bool on_rev = edge_helper::unpack_on_rev(packed_edge);
            bool other_rev = edge_helper::unpack_other_rev(packed_edge);
            bool to_curr = edge_helper::unpack_to_curr(packed_edge);
            if (is_rev != on_rev) {
                other_rev ^= 1;
                to_curr ^= 1;
            }
            if (other_rev == target_rev
                && to_curr == ends_here) {
                edges.remove(i);
                edges.remove(i);
                return true;
            }
        }
    }
    std::cerr << "did NOT REMOVE!?" << std::endl;
    return false;
}

/*
bool node_t::remove_edge_to(const bool& is_rev,
                            const uint64_t& target_id,
                            const bool& target_rev) {
    for (uint64_t i = 0; i < edges.size(); ) {
        uint64_t other_id = edges.at(i++);
        if (other_id
        uint8_t packed_edge = edges.at(i++);
        bool on_rev = edge_helper::unpack_on_rev(packed_edge);
        bool other_rev = edge_helper::unpack_other_rev(packed_edge);
        bool to_curr = edge_helper::unpack_to_curr(packed_edge);
        if (target_rev != on_rev) {
            other_rev ^= 1;
            to_curr ^= 1;
        }
        if (other_id == target_id && other_rev == target_rev) {
            i-=2;
            edges.remove(i);
            edges.remove(i);
            return true;
        }
    }
    return false;
}
*/

void node_t::add_path_step(const uint64_t& path_id, const bool& is_rev,
                           const bool& is_start, const bool& is_end,
                           const uint64_t& prev_id, const uint64_t& prev_rank,
                           const uint64_t& next_id, const uint64_t& next_rank) {
    //std::cerr << "packing " << path_id << " " << is_rev << " " << is_start << " " << is_end << std::endl;
    paths.push_back(path_id);
    paths.push_back(step_type_helper::pack(is_rev, is_start, is_end));
    paths.push_back(prev_id);
    paths.push_back(prev_rank);
    paths.push_back(next_id);
    paths.push_back(next_rank);
}

void node_t::add_path_step(const node_t::step_t& step) {
    add_path_step(step.path_id,
                  step.is_rev,
                  step.is_start,
                  step.is_end,
                  step.prev_id,
                  step.prev_rank,
                  step.next_id,
                  step.next_rank);
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
    uint64_t i = PATH_RECORD_LENGTH*rank;
    uint64_t t = paths.at(i+1);
    return {
        paths.at(i),
        step_type_helper::unpack_is_rev(t),
        step_type_helper::unpack_is_start(t),
        step_type_helper::unpack_is_end(t),
        paths.at(i+2),
        paths.at(i+3),
        paths.at(i+4),
        paths.at(i+5),
    };
}

void node_t::for_each_path_step(const std::function<bool(step_t step)>& func) const {
    uint64_t n_paths = path_count();
    for (uint64_t i = 0; i < n_paths; ++i) {
        func(get_path_step(i));
    }
}

void node_t::set_path_step(const uint64_t& rank, const step_t& step) {
    if (rank >= path_count()) assert(false);
    uint64_t i = PATH_RECORD_LENGTH*rank;
    paths[i] = step.path_id;
    paths[i+1] = step_type_helper::pack(step.is_rev,
                                        step.is_start,
                                        step.is_end);
    paths[i+2] = step.prev_id;
    paths[i+3] = step.prev_rank;
    paths[i+4] = step.next_id;
    paths[i+5] = step.next_rank;
}

bool node_t::step_is_start(const uint64_t& rank) const {
    return step_type_helper::unpack_is_start(paths.at(PATH_RECORD_LENGTH*rank+1));
}

bool node_t::step_is_end(const uint64_t& rank) const {
    return step_type_helper::unpack_is_end(paths.at(PATH_RECORD_LENGTH*rank+1));
}

std::pair<std::map<uint64_t, std::pair<uint64_t, bool>>, // path fronts
          std::map<uint64_t, std::pair<uint64_t, bool>>> // path backs
node_t::flip_paths(void) {
    const std::vector<node_t::step_t> steps = get_path_steps();
    // remove all path steps
    clear_paths();
    // flip them and replace, recording which path starts and ends should be rewritten
    std::pair<std::map<uint64_t, std::pair<uint64_t, bool>>,
              std::map<uint64_t, std::pair<uint64_t, bool>>> path_start_end_rewrites;
    uint64_t rank = 0;
    for (auto& step : steps) {
        // flip the step
        step_t flipped = { step.path_id, !step.is_rev,
                           step.is_start, step.is_end,
                           step.prev_id, step.prev_rank,
                           step.next_id, step.next_rank };
        if (step.is_start) {
            path_start_end_rewrites.first[step.path_id] = std::make_pair(rank, !step.is_rev);
        }
        if (step.is_end) {
            path_start_end_rewrites.second[step.path_id] = std::make_pair(rank, !step.is_rev);
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
        paths.remove(offset);
    }
}

void node_t::clear(void) {
    sequence.clear();
    clear_edges();
    clear_paths();
}

void node_t::clear_edges(void) {
    dyn::hacked_vector null_iv;
    edges = null_iv;
}

void node_t::clear_paths(void) {
    dyn::hacked_vector null_iv;
    paths = null_iv;
}

uint64_t node_t::serialize(std::ostream& out) const {
    uint64_t written = 0;
    size_t seq_size = sequence.size();
    out.write((char*)&seq_size, sizeof(size_t));
    written += sizeof(size_t);
    out.write((char*)sequence.c_str(), seq_size*sizeof(char));
    written += seq_size*sizeof(char);
    written += edges.serialize(out);
    written += paths.serialize(out);
    return written;
}

void node_t::load(std::istream& in) {
    size_t len = 0;
    in.read((char*)&len, sizeof(size_t));
    sequence.resize(len);
    in.read((char*)sequence.c_str(), len*sizeof(uint8_t));
    edges.load(in);
    paths.load(in);
    //display();
}

void node_t::display(void) const {
    std::cerr << "seq " << sequence << " "
              << "edge_count " << edge_count() << " "
              << "path_count " << path_count();
    std::cerr << " | ";
    if (edge_count()) {
        for (uint64_t i = 0; i < edge_count(); ++i) {
            std::cerr
                << edges.at(i) << ":"
                << edges.at(i+1) << " ";
        }
    }
    std::cerr << " | ";
    if (path_count()) {
        for_each_path_step(
            [&](const step_t& step) {
                std::cerr
                    << step.path_id << ":"
                    << step.is_rev << ":"
                    << step.is_start << ":"
                    << step.is_end << ":"
                    << step.prev_id << ":"
                    << step.prev_rank << ":"
                    << step.next_id << ":"
                    << step.next_rank << " ";
                return true;
            });
    }
    std::cerr << std::endl;
}

}
