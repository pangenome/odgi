#pragma once

/**
 * \file gfa_to_handle.hpp
 *
 * Contains a method to construct a mutable handle graph out of a GFA file
 *
 */

#include "gfakluge.hpp"
#include <iostream>
#include <limits>
#include <handlegraph/mutable_path_mutable_handle_graph.hpp>
#include <atomic>
#include <thread>
#include <mutex>
#include <functional>
#include "atomic_queue.h"
#include "progress.hpp"

namespace odgi {

struct path_elem_t {
    handlegraph::path_handle_t path;
    gfak::path_line_t gfak;
};

typedef atomic_queue::AtomicQueue<path_elem_t*, 2 << 10> gfa_path_queue_t;

/// Only the fields gfa_to_handle actually needs from gfak::edge_elem, kept as a
/// plain value type so it can be pushed into gfa_edge_queue_t without a
/// per-edge heap allocation (unlike gfak::edge_elem, which also carries an
/// unused CIGAR alignment string and tags map).
struct edge_record_t {
    std::string source_name;
    std::string sink_name;
    bool source_orientation_forward = false;
    bool sink_orientation_forward = false;
};

typedef atomic_queue::AtomicQueue2<edge_record_t, 2 << 10> gfa_edge_queue_t;

std::map<char, uint64_t> gfa_line_counts(const char* filename);

/// Fills a handle graph with an instantiation of a sequence graph from a GFA file.
/// Handle graph must be empty when passed into function.
void gfa_to_handle(const string& gfa_filename,
                   handlegraph::MutablePathMutableHandleGraph* graph,
                   bool compact_ids,
                   uint64_t n_threads,
                   bool show_progress);

}
