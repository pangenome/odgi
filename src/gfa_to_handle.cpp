#include "gfa_to_handle.hpp"
#if defined(ODGI_PARALLEL_NODE_BUILD) || defined(ODGI_PARALLEL_EDGE_BUILD)
#include "odgi.hpp"
#include <memory>
#endif

namespace odgi {

std::map<char, uint64_t> gfa_line_counts(const char* filename) {
    int gfa_fd = -1;
    char* gfa_buf = nullptr;
    size_t gfa_filesize = gfak::mmap_open(filename, gfa_buf, gfa_fd);
    if (gfa_fd == -1) {
        cerr << "Couldn't open GFA file " << filename << "." << endl;
        exit(1);
    }
    string line;
    size_t i = 0;
    //bool seen_newline = true;
    std::map<char, uint64_t> counts;
    while (i < gfa_filesize) {
        if (i == 0 || gfa_buf[i-1] == '\n') {
            counts[gfa_buf[i]]++;
        }
        ++i;
    }
    gfak::mmap_close(gfa_buf, gfa_fd, gfa_filesize);
    return counts;
}

void gfa_to_handle(const string& gfa_filename,
                   handlegraph::MutablePathMutableHandleGraph* graph,
                   bool compact_ids,
                   uint64_t n_threads,
                   bool progress) {

    n_threads = (n_threads == 0 ? 1 : n_threads);
    char* filename = (char*) gfa_filename.c_str();
    //std::cerr << "filename is " << filename << std::endl;
    gfak::GFAKluge gg;
    uint64_t i = 0;
    uint64_t min_id = std::numeric_limits<uint64_t>::max();
    uint64_t max_id = std::numeric_limits<uint64_t>::min();
    std::map<char, uint64_t> line_counts;

    auto phase_start = std::chrono::steady_clock::now();
    auto log_phase = [&](const char* name) {
        if (progress) {
            auto now = std::chrono::steady_clock::now();
            double secs = std::chrono::duration<double>(now - phase_start).count();
            std::cerr << "[odgi::gfa_to_handle] [timing] " << name << ": " << secs << "s" << std::endl;
            phase_start = now;
        }
    };

    // in parallel scan over the file to count edges and sequences
    {
        std::thread x(
            [&]() {
                gg.for_each_sequence_line_in_file(
                    filename,
                    [&](gfak::sequence_elem s) {
                        try {
                            uint64_t id = stol(s.name);
                            min_id = std::min(min_id, id);
                            max_id = std::max(max_id, id);
                        } catch (const std::exception& e) {
                            std::cerr << "[odgi::gfa_to_handle] Error parsing segment '" << s.name << "': " << e.what() << std::endl;
                            exit(1);
                        }
                    });
            });
        line_counts = gfa_line_counts(filename);
        x.join();
    }
    log_phase("pre-scan (min/max id + line counts)");
    uint64_t id_increment = (compact_ids ? min_id - 1 : 0);
    uint64_t node_count = line_counts['S'];
    uint64_t edge_count = line_counts['L'] + line_counts['E']; // GFA1 'L' and GFA2 'E' edge lines
    uint64_t path_count = line_counts['P'];
    // build the nodes
    {
        std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_meter;
        if (progress) {
            progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                node_count, "[odgi::gfa_to_handle] building nodes:");
        }
        auto build_nodes_serial = [&]() {
            gg.for_each_sequence_line_in_file(
                filename,
                [&](const gfak::sequence_elem& s) {
                    if (s.name.empty() || s.name.find_first_not_of("0123456789") != std::string::npos) {
                        std::cerr << "[odgi::gfa_to_handle] error: segment name '" << s.name
                                  << "' is not a non-negative integer node id" << std::endl;
                        exit(1);
                    }
                    uint64_t id = 0;
                    try {
                        id = std::stoull(s.name);
                    } catch (const std::exception& e) {
                        std::cerr << "[odgi::gfa_to_handle] error: could not parse segment name '" << s.name
                                  << "' as a node id: " << e.what() << std::endl;
                        exit(1);
                    }
                    const uint64_t node_id = id - id_increment;
                    if (graph->has_node(node_id)) {
                        std::cerr << "[odgi::gfa_to_handle] error: duplicate node id " << node_id
                                  << " (segment '" << s.name << "'); GFA node ids must be unique" << std::endl;
                        exit(1);
                    }
                    graph->create_handle(s.sequence, node_id);
                    if (progress) progress_meter->increment(1);
                });
        };
#ifdef ODGI_PARALLEL_NODE_BUILD
        // Experimental (enable via -DODGI_PARALLEL_NODE_BUILD): only takes
        // effect for odgi's own graph_t and n_threads > 1; other handle
        // graph implementations, or n_threads == 1, fall back to
        // build_nodes_serial() above (identical to the always-serial
        // behavior when this flag is off). Requires knowing the id range
        // up front (already available from the pre-scan above) to pre-size
        // storage, so threads can fill distinct, disjoint node_v slots
        // with no shared mutable state -- see graph_t::reserve_node_space
        // and friends in odgi.hpp for the safety argument. Duplicate ids
        // (a hard GFA error) are still detected deterministically via a
        // separate atomic claim array, since reading/writing the same
        // node_v slot from two threads without one would itself be a race.
        odgi::graph_t* fast_graph = node_count > 0 ? dynamic_cast<odgi::graph_t*>(graph) : nullptr;
        if (fast_graph && n_threads > 1) {
            uint64_t max_node_rank = max_id - id_increment;
            fast_graph->reserve_node_space(max_node_rank);
            auto claimed = std::make_unique<std::atomic<bool>[]>(max_node_rank);
            gg.for_each_sequence_endpoints_in_file_parallel(
                filename, n_threads,
                [&](const gfak::sequence_record_t& s) {
                    if (s.name.empty() || s.name.find_first_not_of("0123456789") != std::string::npos) {
                        std::cerr << "[odgi::gfa_to_handle] error: segment name '" << s.name
                                  << "' is not a non-negative integer node id" << std::endl;
                        exit(1);
                    }
                    uint64_t id = 0;
                    try {
                        id = std::stoull(std::string(s.name));
                    } catch (const std::exception& e) {
                        std::cerr << "[odgi::gfa_to_handle] error: could not parse segment name '" << s.name
                                  << "' as a node id: " << e.what() << std::endl;
                        exit(1);
                    }
                    const uint64_t node_id = id - id_increment;
                    if (claimed[node_id - 1].exchange(true, std::memory_order_acq_rel)) {
                        std::cerr << "[odgi::gfa_to_handle] error: duplicate node id " << node_id
                                  << " (segment '" << s.name << "'); GFA node ids must be unique" << std::endl;
                        exit(1);
                    }
                    fast_graph->create_handle_prereserved(std::string(s.sequence), node_id);
                    if (progress) progress_meter->increment(1);
                });
            fast_graph->finalize_prereserved_node_space(min_id - id_increment, max_node_rank);
        } else {
            build_nodes_serial();
        }
#else
        build_nodes_serial();
#endif
        if (progress) {
            progress_meter->finish();
        }
    }
    log_phase("building nodes");

    // building edges and paths: a single parallel pass over the file finds
    // both 'L'/'E' and 'P' lines (gfak::for_each_edge_and_path_line_in_file_parallel),
    // instead of two separate full-file scans. Paths always go through a
    // queue+worker-pool pipeline (append_step needs real graph-structure
    // synchronization). Edges normally do too (edge_queue/edge_worker
    // below) -- except when ODGI_PARALLEL_EDGE_BUILD is enabled and graph
    // is odgi's own graph_t, in which case edges instead go through a
    // lock-free CSR bulk build (see graph_t::bulk_build_edges_from_candidates
    // in odgi.hpp/odgi.cpp for why: the per-edge path here calls
    // create_edge(), which serializes every insertion behind a per-node
    // lock AND does an O(current degree) duplicate-edge scan under that
    // lock -- O(k^2) for a node with k edges, catastrophic for the
    // extreme-degree "hub" nodes real pangenome graphs have at variant-dense
    // loci).
    {
        std::unique_ptr<algorithms::progress_meter::ProgressMeter> edge_progress_meter;
        if (progress) {
            edge_progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                edge_count, "[odgi::gfa_to_handle] building edges:");
        }
        std::unique_ptr<algorithms::progress_meter::ProgressMeter> path_progress_meter;
        if (progress && path_count > 0) {
            path_progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                path_count, "[odgi::gfa_to_handle] building paths:");
        }

        gfa_edge_queue_t edge_queue;
        gfa_path_queue_t path_queue;
        std::atomic<bool> edge_work_todo{};
        std::atomic<bool> path_work_todo{};

        auto edge_worker =
            [&](uint64_t tid) {
                while (edge_work_todo.load()) {
                    edge_record_t e;
                    if (edge_queue.try_pop(e)) {
                        if (e.source_name.empty()) {
                            continue;
                        }
                        try {
                            uint64_t source_id = stol(e.source_name) - id_increment;
                            uint64_t sink_id = stol(e.sink_name) - id_increment;
                            if (graph->has_node(source_id) && graph->has_node(sink_id)) {
                                handlegraph::handle_t a = graph->get_handle(source_id, !e.source_orientation_forward);
                                handlegraph::handle_t b = graph->get_handle(sink_id, !e.sink_orientation_forward);
                                graph->create_edge(a, b);
                            } else {
                                std::cerr << "[odgi::gfa_to_handle] Error creating edge '" << e.source_name << " <--> " << e.sink_name << "' due to missing node(s)" << std::endl;
                                exit(1);
                            }
                        } catch (const std::exception& exc) {
                            std::cerr << "[odgi::gfa_to_handle] Error creating edge '" << e.source_name << " <--> " << e.sink_name << "': " << exc.what() << std::endl;
                            exit(1);
                        }
                        if (progress) edge_progress_meter->increment(1);
                    } else {
                        std::this_thread::sleep_for(std::chrono::nanoseconds(1));
                    }
                }
            };

#ifdef ODGI_PARALLEL_EDGE_BUILD
        // CSR candidate buffer for the fast path. Sized to the worst case
        // (every line touches two distinct node-sides); n_used tracks how
        // many slots were actually claimed (same-rank edges use only one).
        odgi::graph_t* fast_graph_edges = edge_count > 0 ? dynamic_cast<odgi::graph_t*>(graph) : nullptr;
        std::vector<odgi::edge_candidate_t> edge_candidates;
        std::atomic<uint64_t> edge_candidate_next{0};
        if (fast_graph_edges) {
            edge_candidates.resize(edge_count * 2);
        }
        auto collect_edge_candidate =
            [&](const gfak::edge_endpoints_t& e) {
                if (e.source_name.empty()) {
                    return;
                }
                uint64_t source_id, sink_id;
                try {
                    source_id = stol(std::string(e.source_name)) - id_increment;
                    sink_id = stol(std::string(e.sink_name)) - id_increment;
                } catch (const std::exception& exc) {
                    std::cerr << "[odgi::gfa_to_handle] Error creating edge '" << e.source_name << " <--> " << e.sink_name << "': " << exc.what() << std::endl;
                    exit(1);
                }
                if (!(graph->has_node(source_id) && graph->has_node(sink_id))) {
                    std::cerr << "[odgi::gfa_to_handle] Error creating edge '" << e.source_name << " <--> " << e.sink_name << "' due to missing node(s)" << std::endl;
                    exit(1);
                }
                bool source_rev = !e.source_orientation_forward;
                bool sink_rev = !e.sink_orientation_forward;
                // matches graph_t::create_edge's own left/right add_edge calls exactly
                // (see bulk_build_edges_from_candidates doc comment in odgi.hpp)
                uint64_t slot1 = edge_candidate_next.fetch_add(1, std::memory_order_relaxed);
                edge_candidates[slot1] = odgi::edge_candidate_t{source_id - 1, sink_id, sink_rev, false, source_rev};
                if (source_id != sink_id) {
                    uint64_t slot2 = edge_candidate_next.fetch_add(1, std::memory_order_relaxed);
                    edge_candidates[slot2] = odgi::edge_candidate_t{sink_id - 1, source_id, source_rev, true, sink_rev};
                }
                if (progress) edge_progress_meter->increment(1);
            };
#endif

        auto path_worker =
            [&](uint64_t tid) {
                while (path_work_todo.load()) {
                    path_elem_t * p;
                    if (path_queue.try_pop(p)) {
                        uint64_t i = 0;
                        for (auto& s : p->gfak.segment_names) {
                            if (s.empty()) { ++i; continue; } // empty path field: stepless path, don't abort
                            uint64_t id = 0;
                            try {
                                size_t parsed = 0;
                                id = std::stoull(s, &parsed) - id_increment;
                                if (parsed != s.size()) { // reject trailing junk, e.g. a space before * instead of a tab
                                    std::cerr << "[odgi::gfa_to_handle] error: malformed path segment '" << s
                                              << "' in path '" << graph->get_path_name(p->path) << "'" << std::endl;
                                    exit(1);
                                }
                                if (graph->has_node(id)) {
                                    graph->append_step(p->path,
                                                graph->get_handle(id,
                                                                    // in gfak, true == +
                                                                    !p->gfak.orientations[i++]));
                                } else {
                                    std::cerr << "[odgi::gfa_to_handle] Error creating path '" << graph->get_path_name(p->path) << "' due to missing node '" << s << "'" << std::endl;
                                    exit(1);
                                }
                            } catch (...) {
                                std::cerr << "[odgi::gfa_to_handle] id parsing failure for path "
                                          << graph->get_path_name(p->path)
                                          << " attempting to parse node id from '" << s << "'" << std::endl;
                                exit(1);
                            }
                        }
                        delete p;
                        if (progress) path_progress_meter->increment(1);
                    } else {
                        std::this_thread::sleep_for(std::chrono::nanoseconds(1));
                    }
                }
            };

#ifdef ODGI_PARALLEL_EDGE_BUILD
        bool use_csr_edges = (fast_graph_edges != nullptr);
#else
        bool use_csr_edges = false;
#endif

        std::vector<std::thread> edge_workers;
        if (!use_csr_edges) {
            edge_workers.reserve(n_threads);
            edge_work_todo.store(true);
            for (uint64_t t = 0; t < n_threads; ++t) {
                edge_workers.emplace_back(edge_worker, t);
            }
        }

        std::vector<std::thread> path_workers;
        if (path_count > 0) {
            path_workers.reserve(n_threads);
            path_work_todo.store(true);
            for (uint64_t t = 0; t < n_threads; ++t) {
                path_workers.emplace_back(path_worker, t);
            }
        }

        auto path_line_cb =
            [&](const gfak::path_line_t& p) {
                handlegraph::path_handle_t p_h = graph->create_path_handle(p.name);
                path_elem_t* pe = new path_elem_t({p_h, p});
                path_queue.push(pe);
            };

#ifdef ODGI_PARALLEL_EDGE_BUILD
        if (use_csr_edges) {
            gg.for_each_edge_and_path_line_in_file_parallel(
                filename, n_threads, collect_edge_candidate, path_line_cb);
            log_phase("edge+path scan (file read/parse, CSR edge candidates collected inline)");
            if (progress) {
                edge_progress_meter->finish();
            }
            uint64_t n_used = edge_candidate_next.load(std::memory_order_relaxed);
            fast_graph_edges->bulk_build_edges_from_candidates(edge_candidates, n_used, n_threads);
            log_phase("bulk edge build (CSR count/prefix-sum/scatter/finalize)");
        } else
#endif
        {
            gg.for_each_edge_and_path_line_in_file_parallel(
                filename, n_threads,
                [&](const gfak::edge_endpoints_t& e) {
                    edge_queue.push(edge_record_t{std::string(e.source_name), std::string(e.sink_name),
                                                  e.source_orientation_forward, e.sink_orientation_forward});
                },
                path_line_cb);
            log_phase("edge+path scan (file read/parse/enqueue, concurrent with worker drain)");

            while (!edge_queue.was_empty()) {
                std::this_thread::sleep_for(std::chrono::nanoseconds(1));
            }
            edge_work_todo.store(false);
            for (uint64_t t = 0; t < n_threads; ++t) {
                edge_workers[t].join();
            }
            if (progress) {
                edge_progress_meter->finish();
            }
            log_phase("edge worker drain (after scan finished)");
        }

        if (path_count > 0) {
            while (!path_queue.was_empty()) {
                std::this_thread::sleep_for(std::chrono::nanoseconds(1));
            }
            path_work_todo.store(false);
            for (uint64_t t = 0; t < n_threads; ++t) {
                path_workers[t].join();
            }
            if (progress) {
                path_progress_meter->finish();
            }
            log_phase("path worker drain (after scan finished)");
        }
    }

    if (compact_ids) {
        graph->optimize();
        log_phase("optimize (compact_ids)");
    }

}
}
