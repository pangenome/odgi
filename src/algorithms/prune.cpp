#include "prune.hpp"
#include "algorithms/subgraph/extract.hpp" // make_path_name, for consistent subpath naming

namespace odgi {

namespace algorithms {

static_assert(sizeof(recorded_step_t) == 12, "recorded_step_t should pack into 12 bytes");

std::vector<recorded_path_t> record_paths(const graph_t& graph) {
    // The node id is kept full-width; only the length is packed (31 bits), so guard against overflow.
    std::vector<recorded_path_t> recorded;
    recorded.reserve(graph.get_path_count());
    graph.for_each_path_handle([&](const path_handle_t& p) {
        recorded_path_t rp;
        rp.name = graph.get_path_name(p);
        rp.is_circular = graph.get_is_circular(p);
        rp.steps.reserve(graph.get_step_count(p));
        graph.for_each_step_in_path(p, [&](const step_handle_t& step) {
            const handle_t h = graph.get_handle_of_step(step);
            const uint64_t len = graph.get_length(h);
            if (len > 0x7FFFFFFFull) {
                std::cerr << "[odgi::prune] error: node " << graph.get_id(h) << " length " << len
                          << " exceeds 2^31; -S/-A path rebuilding is unsupported for this graph." << std::endl;
                exit(1);
            }
            recorded_step_t s;
            s.set_node_id(graph.get_id(h));
            s.length = (uint32_t)len;
            s.is_rev = graph.get_is_reverse(h);
            rp.steps.push_back(s);
        });
        recorded.push_back(std::move(rp));
    });
    return recorded;
}

prune_path_rebuild_stats_t rebuild_pruned_paths(const std::vector<recorded_path_t>& recorded_paths,
                                                graph_t& subgraph,
                                                const affected_path_policy_t policy) {
    prune_path_rebuild_stats_t stats;
    // A run is a maximal stretch of steps that is a valid walk in the pruned subgraph: every node
    // present, consecutive nodes joined by a surviving edge.
    struct run_t {
        uint64_t start = 0;                 // path offset (bp) at the run's start
        uint64_t end = 0;                   // path offset (bp) at the run's end
        std::vector<handle_t> handles;      // subgraph handles in path order
    };

    // Single-threaded: concurrent append_step to shared nodes races on the in-place graph.
    for (auto& rp : recorded_paths) {
        std::vector<run_t> runs;
        uint64_t walked = 0;      // bp walked along the recorded path
        bool in_run = false;
        handle_t prev_dest;

        for (auto& st : rp.steps) {
            if (subgraph.has_node(st.node_id())) {
                const handle_t dest_h = subgraph.get_handle(st.node_id(), st.is_rev);
                if (!in_run) {
                    runs.push_back({walked, 0, {dest_h}});
                    in_run = true;
                } else if (subgraph.has_edge(prev_dest, dest_h)) {
                    runs.back().handles.push_back(dest_h);
                } else {
                    // both nodes survived but the connecting edge was pruned: break here
                    runs.back().end = walked;
                    runs.push_back({walked, 0, {dest_h}});
                }
                prev_dest = dest_h;
            } else if (in_run) {
                runs.back().end = walked;
                in_run = false;
            }
            walked += st.length;
        }
        if (in_run) {
            runs.back().end = walked;
        }
        const uint64_t path_len = walked;

        // An empty path (no steps) cannot be affected; treat it as intact.
        if (rp.steps.empty()) {
            subgraph.create_path_handle(rp.name, rp.is_circular);
            ++stats.paths_intact;
            continue;
        }

        const bool intact = (runs.size() == 1
                             && runs.front().start == 0
                             && runs.front().end == path_len);

        if (intact) {
            const path_handle_t np = subgraph.create_path_handle(rp.name, rp.is_circular);
            for (auto& h : runs.front().handles) {
                subgraph.append_step(np, h);
            }
            ++stats.paths_intact;
        } else if (policy == affected_path_policy_t::drop_affected) {
            ++stats.paths_dropped;
        } else { // split
            if (runs.empty()) {
                ++stats.paths_dropped;
            } else {
                for (auto& run : runs) {
                    const path_handle_t np = subgraph.create_path_handle(
                        make_path_name(rp.name, run.start, run.end), false);
                    for (auto& h : run.handles) {
                        subgraph.append_step(np, h);
                    }
                    ++stats.subpaths_created;
                }
                ++stats.paths_split;
            }
        }
    }
    return stats;
}

std::vector<edge_t> find_edges_to_prune(const HandleGraph& graph,
                                        size_t k, size_t edge_max,
                                        int n_threads) {
    // for each position on the forward and reverse of the graph
    //unordered_set<edge_t> edges_to_prune;
    std::vector<std::vector<edge_t> > edges_to_prune;
    edges_to_prune.resize(n_threads);
    graph.for_each_handle([&](const handle_t& h) {
            // for the forward and reverse of this handle
            // walk k bases from the end, so that any kmer starting on the node will be represented in the tree we build
            for (auto handle_is_rev : { false, true }) {
                //cerr << "###########################################" << endl;
                handle_t handle = handle_is_rev ? graph.flip(h) : h;
                std::list<walk_t> walks;
                // for each position in the node, set up a kmer with that start position and the node end or kmer length as the end position
                // determine next positions
                nid_t handle_id = graph.get_id(handle);
                size_t handle_length = graph.get_length(handle);
                for (size_t i = 0; i < handle_length;  ++i) {
                    pos_t begin = make_pos_t(handle_id, handle_is_rev, handle_length);
                    pos_t end = make_pos_t(handle_id, handle_is_rev, std::min(handle_length, i+k));
                    walk_t walk = walk_t(offset(end)-offset(begin), begin, end, handle, 0);
                    if (walk.length < k) {
                        // are we branching over more than one edge?
                        size_t next_count = 0;
                        graph.follow_edges(walk.curr, false, [&](const handle_t& next) { ++next_count; return next_count <= 1; });
                        graph.follow_edges(walk.curr, false, [&](const handle_t& next) {
                                if (next_count > 1 && edge_max == walk.forks) { // our next step takes us over the max
                                    int tid = omp_get_thread_num();
                                    edges_to_prune[tid].push_back(graph.edge_handle(walk.curr, next));
                                } else {
                                    walks.push_back(walk);
                                    auto& todo = walks.back();
                                    todo.curr = next;
                                    if (next_count > 1) {
                                        ++todo.forks;
                                    }
                                }
                            });
                    } else {
                        walks.push_back(walk);
                    }
                }
                // now expand the kmers until they reach k
                while (!walks.empty()) {
                    // first we check which ones have reached length k in the current handle; for each of these we run lambda and remove them from our list
                    auto walks_end = walks.end();
                    for (std::list<walk_t>::iterator q = walks.begin(); q != walks_end; ++q) {
                        auto& walk = *q;
                        // did we reach our target length?
                        if (walk.length >= k) {
                            q = walks.erase(q);
                        } else {
                            nid_t curr_id = graph.get_id(walk.curr);
                            size_t curr_length = graph.get_length(walk.curr);
                            bool curr_is_rev = graph.get_is_reverse(walk.curr);
                            size_t take = std::min(curr_length, k-walk.length);
                            walk.end = make_pos_t(curr_id, curr_is_rev, take);
                            walk.length += take;
                            if (walk.length < k) {
                                // if not, we need to expand through the node then follow on
                                size_t next_count = 0;
                                graph.follow_edges(walk.curr, false, [&](const handle_t& next) { ++next_count; return next_count <= 1; });
                                graph.follow_edges(walk.curr, false, [&](const handle_t& next) {
                                        if (next_count > 1 && edge_max == walk.forks) { // our next step takes us over the max
                                            int tid = omp_get_thread_num();
                                            edges_to_prune[tid].push_back(graph.edge_handle(walk.curr, next));
                                        } else {
                                            walks.push_back(walk);
                                            auto& todo = walks.back();
                                            todo.curr = next;
                                            if (next_count > 1) {
                                                ++todo.forks;
                                            }
                                        }
                                    });
                                q = walks.erase(q);
                            } else {
                                // nothing, we'll remove it next time around
                            }
                        }
                    }
                }
            }
        }, true);
    uint64_t total_edges = 0;
    for (auto& v : edges_to_prune) total_edges += v.size();
    std::vector<edge_t> merged; merged.reserve(total_edges);
    for (auto& v : edges_to_prune) {
        merged.insert(merged.end(), v.begin(), v.end());
    }
    // duplicates are assumed to be dealt with externally
    return merged;
}

}

}
