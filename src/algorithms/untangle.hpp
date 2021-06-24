#pragma once

#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <vector>
#include <set>
#include <deque>
#include "hash_map.hpp"

namespace odgi {
namespace algorithms {
    
using namespace handlegraph;

std::vector<step_handle_t> untangle_cuts(
    const PathHandleGraph& graph,
    const step_handle_t& start,
    const step_handle_t& end,
    const ska::flat_hash_map<step_handle_t, uint64_t>& step_pos);

std::vector<step_handle_t> merge_cuts(
    const std::vector<step_handle_t>& cuts,
    const uint64_t& dist,
    const ska::flat_hash_map<step_handle_t, uint64_t>& step_pos);

void write_cuts(
    const PathHandleGraph& graph,
    const path_handle_t& path,
    const std::vector<step_handle_t>& cuts,
    const ska::flat_hash_map<step_handle_t, uint64_t>& path_pos);

ska::flat_hash_map<step_handle_t, uint64_t> make_step_index(
    const PathHandleGraph& graph,
    const std::vector<path_handle_t>& paths);

void untangle(
    const PathHandleGraph& graph,
    const std::vector<path_handle_t>& queries,
    const std::vector<path_handle_t>& targets,
    const uint64_t& merge_dist,
    const size_t& num_threads);

}
}


