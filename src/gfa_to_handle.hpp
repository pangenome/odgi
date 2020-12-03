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

namespace odgi {

struct path_elem_t {
    handlegraph::path_handle_t path;
    gfak::path_elem gfak;
};

typedef atomic_queue::AtomicQueue<path_elem_t*, 2 << 10> gfa_path_queue_t;

/// Fills a handle graph with an instantiation of a sequence graph from a GFA file.
/// Handle graph must be empty when passed into function.
void gfa_to_handle(const string& gfa_filename,
                   handlegraph::MutablePathMutableHandleGraph* graph,
                   uint64_t n_threads = 1,
                   bool show_progress = false);

}
