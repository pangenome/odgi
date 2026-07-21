#include "gfa_graph_to_handle.hpp"
#include "atomic_queue.h"
#include "progress.hpp"

#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace odgi {

// One path (or walk) to append, handed off to a worker thread. `steps` points
// into the decoded GfaGraph, which stays alive until the workers have joined.
struct gfa_graph_path_job_t {
    handlegraph::path_handle_t path;
    const std::vector<NodeId>* steps;
};

typedef atomic_queue::AtomicQueue<gfa_graph_path_job_t*, 2 << 10> gfa_graph_path_queue_t;

// GFAz reconstructs dense, 1-based node IDs, so any reference is valid iff it
// lies in [1, node_count]. This is a cheap range check (no hash lookup) that
// still refuses to silently build a corrupt graph from a malformed .gfaz.
static inline void check_node_ref(uint64_t id, uint64_t node_count, const char* what) {
    if (id < 1 || id > node_count) {
        std::cerr << "[odgi::gfa_graph_to_handle] error: " << what
                  << " references node id " << id
                  << " outside the valid range [1, " << node_count << "]" << std::endl;
        exit(1);
    }
}

// Drain a queue of path/walk jobs, appending each step to its path. Node IDs
// are validated against node_count before use.
static void append_steps_worker(gfa_graph_path_queue_t& queue,
                                std::atomic<bool>& work_todo,
                                uint64_t node_count,
                                handlegraph::MutablePathMutableHandleGraph* graph,
                                algorithms::progress_meter::ProgressMeter* progress_meter) {
    while (work_todo.load()) {
        gfa_graph_path_job_t* job;
        if (queue.try_pop(job)) {
            for (NodeId node_id : *job->steps) {
                uint64_t id = static_cast<uint64_t>(std::abs(node_id));
                check_node_ref(id, node_count, "path/walk step");
                graph->append_step(job->path, graph->get_handle(id, node_id < 0));
            }
            if (progress_meter) progress_meter->increment(1);
            delete job;
        } else {
            std::this_thread::sleep_for(std::chrono::nanoseconds(1));
        }
    }
}

void gfa_graph_to_handle(GfaGraph& gfa_graph,
                         handlegraph::MutablePathMutableHandleGraph* graph,
                         bool compact_ids,
                         uint64_t n_threads,
                         bool show_progress) {

    n_threads = (n_threads == 0 ? 1 : n_threads);

    // node_sequences[0] is a placeholder so that node IDs are 1-based.
    const uint64_t node_count =
        gfa_graph.node_sequences.empty() ? 0 : gfa_graph.node_sequences.size() - 1;

    if (show_progress) {
        std::cerr << "[odgi::gfa_graph_to_handle] building graph from decompressed GFAz: "
                  << node_count << " nodes, "
                  << gfa_graph.links.from_ids.size() << " edges, "
                  << gfa_graph.paths.size() << " paths, "
                  << gfa_graph.walks.size() << " walks (threads=" << n_threads << ")"
                  << std::endl;
    }

    // Columns odgi build never reads. Release them up front to cap peak memory.
    gfa_graph.header_line.clear();
    std::vector<OptionalFieldColumn>().swap(gfa_graph.segment_optional_fields);
    std::vector<OptionalFieldColumn>().swap(gfa_graph.link_optional_fields);
    { JumpData empty; std::swap(gfa_graph.jumps, empty); }
    { ContainmentData empty; std::swap(gfa_graph.containments, empty); }

    // 1. Nodes (S-lines).
    if (node_count > 0) {
        std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_meter;
        if (show_progress) {
            progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                node_count, "[odgi::gfa_graph_to_handle] building nodes:");
        }
        for (uint64_t id = 1; id <= node_count; ++id) {
            // GFAz reconstructs dense 1-based ids, so this cannot trigger in practice; it
            // guards against a malformed .gfaz slipping a duplicate past create_handle's assert.
            if (graph->has_node(id)) {
                std::cerr << "[odgi::gfa_graph_to_handle] error: duplicate node id " << id
                          << " while building nodes from GFAz" << std::endl;
                exit(1);
            }
            graph->create_handle(gfa_graph.node_sequences[id], id);
            if (show_progress) progress_meter->increment(1);
        }
        if (show_progress) progress_meter->finish();
    }

    // Segment names and sequences are done with once the nodes exist.
    std::unordered_map<std::string, uint32_t>().swap(gfa_graph.node_name_to_id);
    std::vector<std::string>().swap(gfa_graph.node_id_to_name);
    std::vector<std::string>().swap(gfa_graph.node_sequences);

    // 2. Edges (L-lines).
    if (!gfa_graph.links.from_ids.empty()) {
        std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_meter;
        if (show_progress) {
            progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                gfa_graph.links.from_ids.size(), "[odgi::gfa_graph_to_handle] building edges:");
        }
        for (size_t i = 0; i < gfa_graph.links.from_ids.size(); ++i) {
            uint64_t source_id = gfa_graph.links.from_ids[i];
            uint64_t sink_id = gfa_graph.links.to_ids[i];
            check_node_ref(source_id, node_count, "edge");
            check_node_ref(sink_id, node_count, "edge");
            handlegraph::handle_t a = graph->get_handle(source_id, gfa_graph.links.from_orients[i] == '-');
            handlegraph::handle_t b = graph->get_handle(sink_id, gfa_graph.links.to_orients[i] == '-');
            graph->create_edge(a, b);
            if (show_progress) progress_meter->increment(1);
        }
        if (show_progress) progress_meter->finish();
    }
    { LinkData empty; std::swap(gfa_graph.links, empty); }

    // 3. Paths (P-lines).
    if (!gfa_graph.paths.empty()) {
        std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_meter;
        if (show_progress) {
            progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                gfa_graph.paths.size(), "[odgi::gfa_graph_to_handle] building paths:");
        }
        gfa_graph_path_queue_t path_queue;
        std::atomic<bool> work_todo{true};
        std::vector<std::thread> workers;
        workers.reserve(n_threads);
        for (uint64_t t = 0; t < n_threads; ++t) {
            workers.emplace_back(append_steps_worker, std::ref(path_queue), std::ref(work_todo),
                                 node_count, graph, progress_meter.get());
        }
        for (size_t i = 0; i < gfa_graph.paths.size(); ++i) {
            handlegraph::path_handle_t p_h = graph->create_path_handle(gfa_graph.path_names[i]);
            path_queue.push(new gfa_graph_path_job_t{p_h, &gfa_graph.paths[i]});
        }
        while (!path_queue.was_empty()) std::this_thread::sleep_for(std::chrono::nanoseconds(1));
        work_todo.store(false);
        for (auto& t : workers) t.join();
        if (show_progress) progress_meter->finish();
    }
    std::vector<std::string>().swap(gfa_graph.path_names);
    std::vector<std::vector<NodeId>>().swap(gfa_graph.paths);
    std::vector<std::string>().swap(gfa_graph.path_overlaps);

    // 4. Walks (W-lines) are intentionally ignored. odgi's text GFA reader
    // (gfa_to_handle) does not import W-lines, so we drop them here too and keep
    // the .gfa and .gfaz builds of the same graph identical. Warn so the dropped
    // data is not silent.
    if (gfa_graph.walks.size() > 0) {
        std::cerr << "[odgi::gfa_graph_to_handle] warning: ignoring "
                  << gfa_graph.walks.size()
                  << " walk(s) (W-lines); odgi does not import walks" << std::endl;
    }
    { WalkData empty; std::swap(gfa_graph.walks, empty); }

    if (compact_ids) {
        graph->optimize();
    }
}

} // namespace odgi
