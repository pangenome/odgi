#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <random>
#include <omp.h>
#include "hash_map.hpp"
#include "position.hpp"
#include "stepindex.hpp"
#include <handlegraph/types.hpp>
#include <handlegraph/iteratee.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <atomic_bitvector.hpp>

namespace odgi {

using namespace handlegraph;

namespace algorithms {

struct step_hash_t {
    step_handle_t step;
    uint64_t pos;
    uint64_t length;
    uint64_t depth;
    uint64_t hash;
    step_hash_t* other_tail = nullptr;
    step_hash_t* other_head = nullptr;
    bool is_head_tail(void) { return other_head != nullptr && other_tail != nullptr; };
};

struct bubble_t {
    step_hash_t* head;
    step_hash_t* tail;
};

struct bubble_chain_t {
    std::vector<bubble_t> bubbles;
};

void for_each_bubble(const PathHandleGraph& graph,
                     const step_index_t& step_index,
                     const path_handle_t& path,
                     const std::function<void(const bubble_t& bubble)>& func);

}

}
