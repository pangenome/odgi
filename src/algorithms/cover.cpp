#include "cover.hpp"

//#define debug_cover

namespace odgi {
    namespace algorithms {

        ska::flat_hash_set<handlegraph::nid_t> is_nice_and_acyclic(const HandleGraph &graph, const ska::flat_hash_set<handlegraph::nid_t> &component) {
            ska::flat_hash_set<handlegraph::nid_t> head_nodes;
            if (component.empty()) { return head_nodes; }

            constexpr size_t NOT_SEEN = std::numeric_limits<size_t>::max();
            std::unordered_map<nid_t, std::pair<size_t, bool>> nodes; // (remaining indegree, orientation)
            std::stack<handle_t> active;
            size_t found = 0; // Number of nodes that have become head nodes.

            // Find the head nodes.
            size_t missing_nodes = 0;
            for (nid_t node : component) {
                if (!(graph.has_node(node))) {
                    missing_nodes++;
                    continue;
                }
                handle_t handle = graph.get_handle(node, false);
                size_t indegree = graph.get_degree(handle, true);
                if (indegree == 0) {
                    nodes[node] = std::make_pair(indegree, false);
                    head_nodes.insert(node);
                    active.push(handle);
                    found++;
                } else {
                    nodes[node] = std::make_pair(NOT_SEEN, false);
                }
            }

            // Active nodes are the current head nodes. Process the successors, determine their
            // orientations, and decrement their indegrees. If the indegree becomes 0, activate
            // the node.
            bool ok = true;
            while (!(active.empty())) {
                handle_t curr = active.top();
                active.pop();
                graph.follow_edges(curr, false, [&](const handle_t &next) -> bool {
                    nid_t next_id = graph.get_id(next);
                    bool next_orientation = graph.get_is_reverse(next);
                    auto iter = nodes.find(next_id);
                    if (iter->second.first == NOT_SEEN) // First visit to the node.
                    {
                        iter->second.first = graph.get_degree(next, true);
                        iter->second.second = next_orientation;
                    } else if (next_orientation != iter->second.second) // Already visited, wrong orientation.
                    {
                        ok = false;
                        return false;
                    }
                    iter->second.first--;
                    if (iter->second.first == 0) {
                        active.push(next);
                        found++;
                    }
                    return true;
                });
                if (!ok) { break; }
            }
            if (found != component.size() - missing_nodes) { ok = false; }

            if (!ok) { head_nodes.clear(); }
            return head_nodes;
        }

        template<class NodeCoverage>
        size_t
        find_first(std::vector<NodeCoverage> &array, nid_t id) {
            size_t first = 0, mid = 0;
            size_t count = array.size();
            while (count > 0) {
                size_t step = count / 2;
                mid = first + step;
                if (array[mid].first < id) {
                    first = mid + 1;
                    count -= step + 1;
                } else { count = step; }
            }
            return first;
        }

        std::vector<handle_t>
        reverse_complement(const HandleGraph &graph, std::vector<handle_t> &forward) {
            std::vector<handle_t> result = forward;
            std::reverse(result.begin(), result.end());
            for (handle_t &handle : result) { handle = graph.flip(handle); }
            return result;
        }

        std::vector<handle_t>
        forward_window(const HandleGraph &graph, const std::deque<handle_t> &path, const handle_t &successor,
                       size_t k) {
            if (path.size() + 1 < k) { k = path.size() + 1; } // Handle the short initial paths in DAGs.
            std::vector<handle_t> forward;
            forward.reserve(k);
            forward.insert(forward.end(), path.end() - (k - 1), path.end());
            forward.push_back(successor);

            std::vector<handle_t> reverse = reverse_complement(graph, forward);
            return (forward < reverse ? forward : reverse);
        }

        std::vector<handle_t>
        backward_window(const HandleGraph &graph, const std::deque<handle_t> &path, const handle_t &predecessor,
                        size_t k) {
            std::vector<handle_t> forward;
            forward.reserve(k);
            forward.push_back(predecessor);
            forward.insert(forward.end(), path.begin(), path.begin() + (k - 1));

            std::vector<handle_t> reverse = reverse_complement(graph, forward);
            return (forward < reverse ? forward : reverse);
        }

        template<class Coverage>
        struct BestCoverage {
            typedef typename Coverage::coverage_t coverage_t;

            coverage_t coverage;
            handle_t handle;

            BestCoverage() : coverage(Coverage::worst_coverage()), handle() {}

            void update(const coverage_t &new_coverage, const handle_t &new_handle) {
                if (new_coverage < this->coverage) {
                    this->coverage = new_coverage;
                    this->handle = new_handle;
                }
            }
        };

        struct SimpleCoverage {
            typedef HandleGraph graph_t;

            typedef size_t coverage_t;
            typedef std::pair<nid_t, coverage_t> node_coverage_t;

            static std::vector<node_coverage_t>
            init_node_coverage(const MutablePathDeletableHandleGraph &graph, const ska::flat_hash_set<handlegraph::nid_t> &component, const bool& ignore_paths) {
                std::vector<node_coverage_t> node_coverage;
                node_coverage.reserve(component.size());
                for (nid_t id : component) {
                    // We are actually interested in the intersection of this graph and the component.
                    // For example, some nodes of the original graph may be missing from a GBWTGraph.
                    if (!(graph.has_node(id))) { continue; }

                    uint64_t depth = ignore_paths ? 0 : graph.get_step_count(graph.get_handle(id));
                    node_coverage.emplace_back(id, static_cast<coverage_t>(depth));
                }

                return node_coverage;
            }

            static bool extend_forward(const graph_t &graph, std::deque<handle_t> &path, size_t k,
                                       std::vector<node_coverage_t> &node_coverage,
                                       std::map<std::vector<handle_t>, coverage_t> &path_coverage, bool acyclic) {
                bool success = false;
                BestCoverage<SimpleCoverage> best;
                graph.follow_edges(path.back(), false, [&](const handle_t &next) {
                    success = true;
                    if (!acyclic && path.size() + 1 < k) // Node coverage.
                    {
                        size_t first = find_first(node_coverage, graph.get_id(next));
                        best.update(node_coverage[first].second, next);
                    } else {
                        std::vector<handle_t> window = forward_window(graph, path, next, k);
                        best.update(path_coverage[window], next);
                    }
                });

                if (success) {
                    if (acyclic || path.size() + 1 >= k) {
                        std::vector<handle_t> window = forward_window(graph, path, best.handle, k);
                        increase_coverage(path_coverage[window]);
                    }
                    if (!acyclic) {
                        size_t first = find_first(node_coverage, graph.get_id(best.handle));
                        increase_coverage(node_coverage[first]);
                    }
                    path.push_back(best.handle);
                }

                return success;
            }

            static bool extend_backward(const graph_t &graph, std::deque<handle_t> &path, size_t k,
                                        std::vector<node_coverage_t> &node_coverage,
                                        std::map<std::vector<handle_t>, coverage_t> &path_coverage) {
                bool success = false;
                BestCoverage<SimpleCoverage> best;
                graph.follow_edges(path.front(), true, [&](const handle_t &prev) {
                    success = true;
                    if (path.size() + 1 < k) // Node coverage.
                    {
                        size_t first = find_first(node_coverage, graph.get_id(prev));
                        best.update(node_coverage[first].second, prev);
                    } else {
                        std::vector<handle_t> window = backward_window(graph, path, prev, k);
                        best.update(path_coverage[window], prev);
                    }
                });

                if (success) {
                    if (path.size() + 1 >= k) {
                        std::vector<handle_t> window = backward_window(graph, path, best.handle, k);
                        increase_coverage(path_coverage[window]);
                    }
                    size_t first = find_first(node_coverage, graph.get_id(best.handle));
                    increase_coverage(node_coverage[first]);
                    path.push_front(best.handle);
                }

                return success;
            }

            static void increase_coverage(coverage_t &coverage) {
                coverage++;
            }

            static void increase_coverage(node_coverage_t &node) {
                increase_coverage(node.second);
            }

            static coverage_t worst_coverage() { return std::numeric_limits<coverage_t>::max(); }

            static std::string name() { return "SimpleCoverage"; }
        };

        template<class Coverage>
        bool
        component_path_cover(handlegraph::MutablePathDeletableHandleGraph &graph,
                             std::vector<ska::flat_hash_set<handlegraph::nid_t>> &components, size_t component_id,
                             size_t num_paths_per_component, size_t node_window_size,
                             size_t min_node_depth, size_t max_number_of_paths_generable,
                             bool write_node_depth, std::string &node_depth,
                             const uint64_t& nthreads, const bool& ignore_paths, const bool& show_progress) {
            typedef typename Coverage::coverage_t coverage_t;
            typedef typename Coverage::node_coverage_t node_coverage_t;

            ska::flat_hash_set<handlegraph::nid_t> &component = components[component_id];
            size_t component_size = component.size();
            ska::flat_hash_set<handlegraph::nid_t> head_nodes = is_nice_and_acyclic(graph, component);
            bool acyclic = !(head_nodes.empty());
            if (show_progress) {
                std::cerr << Coverage::name() << ": processing component " << (component_id + 1) << " / "
                          << components.size() << ", which has " << component.size() << " node(s) and it is "
                          << (acyclic ? "acyclic." : "cyclic.") << std::endl;
            }

            // Node depth for the potential starting nodes.
            std::vector<node_coverage_t> all_nodes_depth = Coverage::init_node_coverage(
                    graph, (acyclic ? head_nodes : component), ignore_paths
                    );
            std::map<std::vector<handle_t>, coverage_t> path_coverage; // Path and its reverse complement are equivalent.

            // all_nodes_depth will be empty if we cannot create this type of path cover for the component.
            // For example, if there are no haplotypes for LocalHaplotypes.
            if (all_nodes_depth.empty()) {
                if (show_progress) {
                    std::cerr << Coverage::name() << ": Cannot find this type of path cover for the component"
                              << std::endl;
                }
                return false;
            }

            // Now that we know that we can find a path cover, we can save a little bit of memory by deleting
            // the component.
            component = ska::flat_hash_set<handlegraph::nid_t>();

            if (num_paths_per_component > 0) {
                min_node_depth = std::numeric_limits<uint64_t>::max();
            } else {
                num_paths_per_component = max_number_of_paths_generable;
            }
            // Generate num_paths_per_component paths in the component.
            uint64_t i;
            for (i = 0; i < num_paths_per_component; i++) {
                // Choose a starting node with minimum coverage and then sort the nodes by id.
                std::deque<handle_t> path;
                ips4o::parallel::sort(all_nodes_depth.begin(), all_nodes_depth.end(),
                          [](const node_coverage_t &a, const node_coverage_t &b) -> bool {
                              return (a.second < b.second);
                          }, nthreads);

#ifdef debug_cover
                std::cerr << node_coverage.front().first << " --- " << node_coverage.front().second << std::endl;
#endif

                if (all_nodes_depth.front().second >= min_node_depth) {
                    if (show_progress) {
                        std::cerr << Coverage::name() << ": minimum node depth reached after generating " << i << " paths." << std::endl;
                    }

                    break;
                }

                path.push_back(graph.get_handle(all_nodes_depth.front().first, false));
                Coverage::increase_coverage(all_nodes_depth.front());
                ips4o::parallel::sort(all_nodes_depth.begin(), all_nodes_depth.end(),
                          [](const node_coverage_t &a, const node_coverage_t &b) -> bool {
                              return (a.first < b.first);
                          }, nthreads);

                // Extend the path forward if acyclic or in both directions otherwise.
                bool success = true;
                while (success && path.size() < component_size) {
                    success = false;
                    success |= Coverage::extend_forward(graph, path, node_window_size, all_nodes_depth, path_coverage,
                                                        acyclic);
                    if (!acyclic && path.size() < component_size) {
                        success |= Coverage::extend_backward(graph, path, node_window_size, all_nodes_depth,
                                                             path_coverage);
                    }
                }

                path_handle_t new_path = graph.create_path_handle(
                        "Path_" + std::to_string(component_id) + "_" + std::to_string(i));
                for (handle_t handle : path) {
                    graph.append_step(new_path, handle);
                }

#ifdef debug_cover
                std::cerr << "Path_" + std::to_string(component_id) + "_" + std::to_string(i) << ":";
                for (handle_t handle : path) {
                    std::cerr << " " << graph.get_id(handle) << (number_bool_packing::unpack_bit(handle) ? "-" : "+");
                }
                std::cerr << std::endl;
#endif
            }

            if ((min_node_depth != std::numeric_limits<uint64_t>::max()) && (i >= num_paths_per_component)){
                std::cerr << Coverage::name() <<": maximum number of generable paths reached." << std::endl;
            }

            if (write_node_depth) {
                if (component_id == 0) {
                    node_depth += "component_id\tnode_id\tdepth\n";
                }

                for (node_coverage_t single_depth : all_nodes_depth) {
                    node_depth +=
                            std::to_string(component_id) + "\t" + std::to_string(single_depth.first) + "\t" +
                            std::to_string(single_depth.second) + "\n";
                }
            }

            return true;
        }

        void path_cover(handlegraph::MutablePathDeletableHandleGraph &graph,
                        size_t num_paths_per_component, size_t node_window_size,
                        size_t min_node_depth, size_t max_number_of_paths_generable,
                        bool write_node_depth, std::string &node_depth,
                        const uint64_t& nthreads, const bool& ignore_paths, const bool& show_progress) {
            std::vector<ska::flat_hash_set<handlegraph::nid_t>> weak_components = algorithms::weakly_connected_components(
                    &graph);

            // Handle each component separately.
            size_t processed_components = 0;
            for (size_t contig = 0; contig < weak_components.size(); contig++) {
                if (component_path_cover<SimpleCoverage>(graph, weak_components, contig,
                                                         num_paths_per_component, node_window_size,
                                                         min_node_depth, max_number_of_paths_generable,
                                                         write_node_depth, node_depth,
                                                         nthreads, ignore_paths, show_progress)) {
                    processed_components++;

                    if (show_progress) {
                        std::cerr << "[odgi::path_cover] Processed: " << processed_components << std::endl;
                    }
                }
            }
        }



void hogwild_path_cover(handlegraph::MutablePathDeletableHandleGraph &graph,
                        double target_depth,
                        const uint64_t& nthreads, const bool& ignore_paths, const bool& show_progress) {

    // get depth -> graph depth
    uint64_t node_count = graph.get_node_count();
    uint64_t graph_bp = 0;
    uint64_t step_count = 0;
    // make a rank select dictionary over our sequence space
    // to get a node randomly distributed in the total length of the graph
    graph.for_each_handle(
        [&](const handle_t& h) {
            graph_bp += graph.get_length(h);
            step_count += graph.get_step_count(h);
        });
    //sdsl::bit_vector graph_vec();
    sdsl::bit_vector graph_bv(graph_bp);
    graph_bp = 0;
    graph.for_each_handle(
        [&](const handle_t& h) {
            graph_bv[graph_bp] = 1;
            graph_bp += graph.get_length(h);
        });
    sdsl::bit_vector::rank_1_type graph_bv_rank;
    sdsl::util::assign(graph_bv_rank, sdsl::bit_vector::rank_1_type(&graph_bv));
    uint64_t target_step_count = node_count * target_depth - (ignore_paths ? 0 : step_count);
    // run hogwild until our thread count is
    std::unique_ptr<progress_meter::ProgressMeter> progress_meter;
    if (show_progress) {
        progress_meter = std::make_unique<progress_meter::ProgressMeter>(
            target_step_count, "[odgi::hogwild_cover] covering the graph:");
    }
    // launch a thread to update the learning rate, count iterations, and decide when to stop
    std::atomic<uint64_t> added_steps; added_steps.store(0);
    std::atomic<uint64_t> added_paths; added_paths.store(0);
    auto worker_lambda =
        [&](uint64_t tid) {
            // everyone tries to seed with their own random data
            const std::uint64_t seed = 9399220 + tid;
            XoshiroCpp::Xoshiro256Plus gen(seed); // a nice, fast PRNG

            // we'll sample from all path steps
            std::uniform_int_distribution<uint64_t> dis_graph_pos = std::uniform_int_distribution<uint64_t>(0, graph_bp-1);
            std::uniform_int_distribution<uint64_t> flip(0, 1);
            while (added_steps.load() < target_step_count) {
                // create a path
                std::stringstream ss;
                ss << "cover_" << added_paths++;
                path_handle_t path = graph.create_path_handle(ss.str());
                // find a random handle
                handle_t h = graph.get_handle(graph_bv_rank(dis_graph_pos(gen)), flip(gen));
                uint64_t iter = 0;
                handle_t lowest;
                bool seen_low = false;
                while (graph.get_step_count(h) > 0 && iter++ < 100) {
                    h = graph.get_handle(graph_bv_rank(dis_graph_pos(gen)), flip(gen));
                    if (!seen_low || graph.get_step_count(h) < graph.get_step_count(lowest)) {
                        lowest = h;
                        seen_low = true;
                    }
                }
                if (seen_low) {
                    h = lowest;
                }
                while (true) {
                    graph.append_step(path, h);
                    if (show_progress) progress_meter->increment(1);
                    ++added_steps;
                    if (added_steps.load() >= target_step_count) {
                        break;
                    }
                    handle_t best_next;
                    uint64_t lowest_cov = std::numeric_limits<uint64_t>::max();
                    graph.follow_edges(
                        h, false,
                        [&](const handle_t& n) {
                            uint64_t next_cov = graph.get_step_count(n);
                            if (next_cov < lowest_cov) {
                                best_next = n;
                                lowest_cov = next_cov;
                            }
                        });
                    if (lowest_cov < std::numeric_limits<uint64_t>::max()) {
                        h = best_next;
                    } else {
                        break;
                    }
                }
            }
        };

    std::vector<std::thread> workers;
    workers.reserve(nthreads);
    for (uint64_t t = 0; t < nthreads; ++t) {
        workers.emplace_back(worker_lambda, t);
    }

    for (uint64_t t = 0; t < nthreads; ++t) {
        workers[t].join();
    }

    if (show_progress) {
        progress_meter->finish();
    }
}


    }
}
