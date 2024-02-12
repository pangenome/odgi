#ifndef DSGVG_POSITION_HPP_INCLUDED
#define DSGVG_POSITION_HPP_INCLUDED

//#include "utility.hpp"
#include <iostream>
#include <tuple>
#include <handlegraph/util.hpp>

/** \file
 * Functions for working with Positions and `pos_t`s.
 */

namespace odgi {

//using namespace handlegraph;

// Resolve ambiguous nid_t typedef by putting it in our namespace.
using nid_t = handlegraph::nid_t;

/// Represents an oriented position on a Node.
/// Position type: id, direction, offset.
/// Offset is counted as for as prorobuf Position, from the node's first base
/// on the forward strand, and from its last base on the reverse strand.
typedef std::tuple<nid_t, bool, uint64_t> pos_t;

/// Make a position that refers to a node id an offset on the node
inline pos_t make_pos_t(nid_t id, bool is_rev, uint64_t off) {
    return std::make_tuple(id, is_rev, off);
}

/// Extract the id of the node a pos_t is on.
inline nid_t id(const pos_t& pos) {
    return std::get<0>(pos);
}

/// Return true if a pos_t is on the reverse strand of its node.
inline bool is_rev(const pos_t& pos) {
    return std::get<1>(pos);
}

/// Get the offset along the selected strand of the node from a pos_t.
inline uint64_t offset(const pos_t& pos) {
    return std::get<2>(pos);
}

/// Get a reference to the Node ID of a pos_t.
inline nid_t& get_id(pos_t& pos) {
    return std::get<0>(pos);
}

/// Get a reference to the reverse flag of a pos_t.
inline bool& get_is_rev(pos_t& pos) {
    return std::get<1>(pos);
}

/// Get a reference to the offset field of a pos_t, which counts along the selected strand of the node.
inline uint64_t& get_offset(pos_t& pos) {
    return std::get<2>(pos);
}

/// Return true if a pos_t is unset.
inline bool is_empty(const pos_t& pos) {
    return id(pos) == 0;
}

/// Reverse a pos_t and get a pos_t at the same base, going the other direction.
inline pos_t reverse(const pos_t& pos, size_t node_length) {
    pos_t rev = pos;
    // swap the offset onto the other strand
    get_offset(rev) = node_length - offset(rev);
    // invert the position
    get_is_rev(rev) = !is_rev(rev);
    return rev;
}

/// A path position
struct path_pos_t {
    handlegraph::path_handle_t path;
    uint64_t offset;
    bool is_rev;
};

struct path_range_t {
    path_pos_t begin;
    path_pos_t end;
    bool is_rev;
    std::string name;
    std::string data;
};

struct path_range_comparator {
    bool operator() (const path_range_t& lhs, const path_range_t& rhs) const {
        if (lhs.begin.path != rhs.begin.path) return lhs.begin.path < rhs.begin.path;
        if (lhs.end.path != rhs.end.path) return lhs.end.path < rhs.end.path;
        if (lhs.begin.offset != rhs.begin.offset) return lhs.begin.offset < rhs.begin.offset;
        return lhs.end.offset < rhs.end.offset;
    }
};

inline std::string& get_long_path_name(std::tuple<std::string, uint64_t, uint64_t> path_long_start_end) {
	return std::get<0>(path_long_start_end);
}

inline uint64_t & get_long_path_start(std::tuple<std::string, uint64_t, uint64_t> path_long_start_end) {
	return std::get<1>(path_long_start_end);
}

inline uint64_t & get_long_path_end(std::tuple<std::string, uint64_t, uint64_t> path_long_start_end) {
	return std::get<2>(path_long_start_end);
}

typedef std::pair<uint64_t, uint64_t> interval_t;

}

#endif
