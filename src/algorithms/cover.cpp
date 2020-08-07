#include "cover.hpp"

//#define debug_cover

namespace odgi {
    namespace algorithms {

        //TODO: check if there is already an implementation like that somewhere
        ska::flat_hash_set<handlegraph::nid_t>
        is_nice_and_acyclic(const HandleGraph &graph, const ska::flat_hash_set<handlegraph::nid_t> &component) {
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
            init_node_coverage(const graph_t &graph, const ska::flat_hash_set<handlegraph::nid_t> &component) {
                std::vector<node_coverage_t> node_coverage;
                node_coverage.reserve(component.size());
                for (nid_t id : component) {
                    // We are actually interested in the intersection of this graph and the component.
                    // For example, some nodes of the original graph may be missing from a GBWTGraph.
                    if (!(graph.has_node(id))) { continue; }
                    node_coverage.emplace_back(id, static_cast<coverage_t>(0));
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
        component_path_cover(handlegraph::MutablePathDeletableHandleGraph& graph,
                             std::vector<ska::flat_hash_set<handlegraph::nid_t>> &components, size_t component_id,
                             size_t n, size_t k, bool show_progress, size_t sample_id_offset, size_t contig_id) {
            typedef typename Coverage::coverage_t coverage_t;
            typedef typename Coverage::node_coverage_t node_coverage_t;

            ska::flat_hash_set<handlegraph::nid_t> &component = components[component_id];
            size_t component_size = component.size();
            ska::flat_hash_set<handlegraph::nid_t> head_nodes = is_nice_and_acyclic(graph, component);
            bool acyclic = !(head_nodes.empty());
            if (show_progress) {
                std::cerr << Coverage::name() << ": Processing component " << (component_id + 1) << " / "
                          << components.size();
                if (acyclic) { std::cerr << " (acyclic)"; }
                std::cerr << std::endl;
            }

            // Node coverage for the potential starting nodes.
            std::vector<node_coverage_t> node_coverage = Coverage::init_node_coverage(graph, (acyclic ? head_nodes
                                                                                                      : component));
            std::map<std::vector<handle_t>, coverage_t> path_coverage; // Path and its reverse complement are equivalent.

            // Node coverage will be empty if we cannot create this type of path cover for the component.
            // For example, if there are no haplotypes for LocalHaplotypes.
            if (node_coverage.empty()) {
                if (show_progress) {
                    std::cerr << Coverage::name() << ": Cannot find this type of path cover for the component"
                              << std::endl;
                }
                return false;
            }

            // Now that we know that we can find a path cover, we can save a little bit of memory by deleting
            // the component.
            component = ska::flat_hash_set<handlegraph::nid_t>();

            // Generate n paths in the component.
            for (size_t i = 0; i < n; i++) {
                // Choose a starting node with minimum coverage and then sort the nodes by id.
                std::deque<handle_t> path;
                std::sort(node_coverage.begin(), node_coverage.end(),
                          [](const node_coverage_t &a, const node_coverage_t &b) -> bool {
                              return (a.second < b.second);
                          });
                path.push_back(graph.get_handle(node_coverage.front().first, false));
                Coverage::increase_coverage(node_coverage.front());
                std::sort(node_coverage.begin(), node_coverage.end(),
                          [](const node_coverage_t &a, const node_coverage_t &b) -> bool {
                              return (a.first < b.first);
                          });

                // Extend the path forward if acyclic or in both directions otherwise.
                bool success = true;
                while (success && path.size() < component_size) {
                    success = false;
                    success |= Coverage::extend_forward(graph, path, k, node_coverage, path_coverage, acyclic);
                    if (!acyclic && path.size() < component_size) {
                        success |= Coverage::extend_backward(graph, path, k, node_coverage, path_coverage);
                    }
                }

                path_handle_t new_path = graph.create_path_handle("Path_" + std::to_string(component_id) + "_" + std::to_string(i));
                for (handle_t handle : path) {
                    graph.append_step(new_path, handle);
                }

#ifdef debug_cover
                std::cerr << "Path_" + std::to_string(component_id) + "_" + std::to_string(i) << " -->";
                for (handle_t handle : path) {
                    std::cerr << " " << graph.get_id(handle) << (number_bool_packing::unpack_bit(handle) ? "-" : "+");
                }
                std::cerr << std::endl;
#endif
            }

            return true;
        }

        void path_cover(handlegraph::MutablePathDeletableHandleGraph &graph,
                        size_t n, size_t k,
                        bool show_progress) {
            std::vector<ska::flat_hash_set<handlegraph::nid_t>> weak_components = algorithms::weakly_connected_components(
                    &graph);

            // Handle each component separately.
            size_t processed_components = 0;
            for (size_t contig = 0; contig < weak_components.size(); contig++) {
                auto &weak_component = weak_components[contig];

#ifdef debug_cover
                std::cerr << "component " << contig << "; size " << weak_component.size() << std::endl;
#endif

                if (component_path_cover<SimpleCoverage>(graph, weak_components, contig, n, k, show_progress, 0,
                                                         contig)) {
                    processed_components++;

                    std::cerr << "Processed: " << processed_components << std::endl;
                }
            }
        }

    }
}
