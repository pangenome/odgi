#include "stepindex.hpp"
#include "progress.hpp"

namespace odgi {
namespace algorithms {

step_index_t::step_index_t(const PathHandleGraph& graph,
                           const std::vector<path_handle_t>& paths,
                           const uint64_t& nthreads,
                           const bool progress) {
    // iterate through the paths, recording steps in the structure we'll use to build the mphf
    std::vector<step_handle_t> steps;
#pragma omp parallel for schedule(dynamic,1)
    for (auto& path : paths) {
        std::vector<step_handle_t> my_steps;
        graph.for_each_step_in_path(
            path, [&](const step_handle_t& step) { my_steps.push_back(step); });
        my_steps.push_back(graph.path_end(path));
#pragma omp critical (steps_collect)
        steps.insert(steps.end(), my_steps.begin(), my_steps.end());
    }
    // sort the steps
    ips4o::parallel::sort(steps.begin(), steps.end(), std::less<>(), nthreads);
    // build the hash function (quietly)
    step_mphf = new boophf_step_t(steps.size(), steps, nthreads, 2.0, false, false);
    // use the hash function to record the step positions
    pos.resize(steps.size());
	std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_meter;
	if (progress) {
		progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
				paths.size(), "[odgi::algorithms::stepindex] Building Progress:");
	}
#pragma omp parallel for schedule(dynamic,1)
    for (auto& path : paths) {
        uint64_t offset = 0;
        graph.for_each_step_in_path(
            path, [&](const step_handle_t& step) {
                pos[step_mphf->lookup(step)] = offset;
                offset += graph.get_length(graph.get_handle_of_step(step));
            });
        pos[step_mphf->lookup(graph.path_end(path))] = offset;
        if (progress) {
        	progress_meter->increment(1);
        }
    }
	if (progress) {
		progress_meter->finish();
	}
}

const uint64_t& step_index_t::get_position(const step_handle_t& step) const {
    return pos[step_mphf->lookup(step)];
}

step_index_t::~step_index_t(void) {
    delete step_mphf;
}


// path step index

path_step_index_t::path_step_index_t(const PathHandleGraph& graph,
                                     const path_handle_t& path,
                                     const uint64_t& nthreads) {
    // iterate through the paths, recording steps in the structure we'll use to build the mphf
    {
        std::vector<step_handle_t> steps;
        std::vector<nid_t> nodes;
        graph.for_each_step_in_path(
            path, [&](const step_handle_t& step) {
                steps.push_back(step);
                nodes.push_back(graph.get_id(graph.get_handle_of_step(step)));
            });
        steps.push_back(graph.path_end(path));
        // sort the steps, nb. they're unique
        ips4o::parallel::sort(steps.begin(), steps.end(), std::less<>(), nthreads);
        // build the hash function (quietly)
        step_mphf = new boophf_step_t(steps.size(), steps, nthreads, 2.0, false, false);

        ips4o::parallel::sort(nodes.begin(), nodes.end(), std::less<>(), nthreads);
        // then take unique positions
        nodes.erase(std::unique(nodes.begin(),
                                nodes.end()),
                    nodes.end());
        // build the hash function (quietly)
        node_mphf = new boophf_uint64_t(nodes.size(), nodes, nthreads, 2.0, false, false);
        node_count = nodes.size();
        step_count = steps.size();
    }

    {
        // here, we sort steps by handle and then position
        // and build our handle->step list and step->offset maps
        // these are steps sorted by the bbhash of their node id, and then their offset in the path
        std::vector<std::tuple<uint64_t, uint64_t, step_handle_t>> steps_by_node;
        uint64_t offset = 0;
        graph.for_each_step_in_path(
            path, [&](const step_handle_t& step) {
                steps_by_node.push_back(std::make_tuple(node_mphf->lookup(graph.get_id(graph.get_handle_of_step(step))),
                                                        offset,
                                                        step));
                offset += graph.get_length(graph.get_handle_of_step(step));
            });

        node_offset.resize(node_count+1);
        step_offset.resize(step_count+1);
        ips4o::parallel::sort(steps_by_node.begin(), steps_by_node.end(), std::less<>(), nthreads);
        uint64_t last_idx = 0;
        node_offset[0] = 0; // first offset is 0 by definition
        for (auto& node_step : steps_by_node) {
            auto& idx = std::get<0>(node_step);
            //auto& offset = std::get<1>(node_step); // just used for sorting
            auto& step = std::get<2>(node_step);
            //std::cerr << "idx = " << idx << " " << as_integers(step)[0] << ":" << as_integers(step)[1] << std::endl;
            if (idx != last_idx) {
                node_offset[idx] = node_steps.size();
            }
            step_offset[step_mphf->lookup(step)] = node_steps.size();
            node_steps.push_back(step);
            last_idx = idx;
        }
        if (last_idx != node_count-1) {
            std::cerr << "last_idx vs node count " << last_idx << " " << node_count << std::endl;
            std::cerr << "[odgi::algorithms::stepindex] unexpected mismatch between last_idx and node_count" << std::endl;
            std::abort();
        }
        node_offset[node_count] = node_steps.size();
        step_offset[step_count] = node_steps.size();
    }
}

path_step_index_t::~path_step_index_t(void) {
    delete node_mphf;
    delete step_mphf;
}

uint64_t path_step_index_t::get_node_idx(const nid_t& id) const {
    return node_mphf->lookup(id);
}

uint64_t path_step_index_t::get_step_idx(const step_handle_t& step) const {
    return step_mphf->lookup(step);
}

uint64_t path_step_index_t::n_steps_on_node(const nid_t& id) const {
    auto idx = get_node_idx(id);
    return node_offset[idx+1] - node_offset[idx];
}

std::pair<bool, step_handle_t>
path_step_index_t::get_next_step_on_node(const nid_t& id, const step_handle_t& step) const {
    auto node_idx = get_node_idx(id);
    auto curr_steps = node_offset[node_idx];
    auto next_steps = node_offset[node_idx+1];
    auto step_idx = step_offset[get_step_idx(step)];
    bool has_next = step_idx + 1 < next_steps;
    if (has_next) {
        return std::make_pair(true, node_steps[step_idx+1]);
    } else {
        step_handle_t empty_step;
        return std::make_pair(false, empty_step);
    }
}

std::pair<bool, step_handle_t>
path_step_index_t::get_prev_step_on_node(const nid_t& id, const step_handle_t& step) const {
    auto curr_steps = node_offset[get_node_idx(id)];
    auto step_idx = step_offset[get_step_idx(step)];
    bool has_prev = step_idx > curr_steps;
    if (has_prev) {
        return std::make_pair(true, node_steps[step_idx-1]);
    } else {
        step_handle_t empty_step;
        return std::make_pair(false, empty_step);
    }
}

}
}
