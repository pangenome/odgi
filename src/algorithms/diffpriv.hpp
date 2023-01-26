#pragma once

#include <vector>
#include <set>
#include <deque>
#include <random>
#include <iostream>
#include <map>
#include <atomic>
#include <thread>
#include <sstream>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
#include <handlegraph/util.hpp>
#include "progress.hpp"
//#include "hash_map.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

void diff_priv(
    const PathHandleGraph& graph,
    MutablePathDeletableHandleGraph& priv,
    //PathHandleGraph& priv,
    const double epsilon,
    const double target_coverage,
    const double min_haplotype_freq,
    const uint64_t bp_limit,
    const uint64_t nthreads,
    const bool progress_reporting,
    const bool write_samples);

// for selecting one-of-a-set-of-steps, adapted from https://gist.github.com/cbsmith/5538174
template <typename RandomGenerator = std::default_random_engine>
struct random_selector {
    random_selector(RandomGenerator g = RandomGenerator(std::random_device()())) : gen(g) {}
    template <typename Iter>
    Iter select(Iter start, Iter end) {
        std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
        std::advance(start, dis(gen));
        return start;
    }
	template <typename Iter>
	Iter operator()(Iter start, Iter end) {
		return select(start, end);
	}
    template <typename Container>
    auto operator()(const Container& c) -> decltype(*begin(c))& {
        return *select(begin(c), end(c));
    }
private:
    RandomGenerator gen;
};

}
}
