#include "topological_sort.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

std::vector<handle_t> head_nodes(const HandleGraph* g) {
    std::vector<handle_t> to_return;
    g->for_each_handle([&](const handle_t& found) {
        // For each (locally forward) node
        
        bool no_left_edges = true;
        g->follow_edges(found, true, [&](const handle_t& ignored) {
            // We found a left edge!
            no_left_edges = false;
            // We only need one
            return false;
        });
        
        if (no_left_edges) {
            to_return.push_back(found);
        }
    });
    
    return to_return;
    
}

std::vector<handle_t> tail_nodes(const HandleGraph* g) {
    std::vector<handle_t> to_return;
    g->for_each_handle([&](const handle_t& found) {
        // For each (locally forward) node
        
        bool no_right_edges = true;
        g->follow_edges(found, false, [&](const handle_t& ignored) {
            // We found a right edge!
            no_right_edges = false;
            // We only need one
            return false;
        });
        
        if (no_right_edges) {
            to_return.push_back(found);
        }
    });
    
    return to_return;
    
}

std::vector<handle_t> topological_order(const HandleGraph* g, bool use_heads, bool use_tails, bool progress_reporting) {

    // Make a vector to hold the ordered and oriented nodes.
    std::vector<handle_t> sorted;
    sorted.reserve(g->get_node_count());
    
    // Instead of actually removing edges, we add them to this set of masked edges.
    dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,256,16> > masked_edges_bv;
    dyn::lciv<dyn::hacked_vector,256,16> masked_edges_iv;

    auto delta_to_handle = [&](const handle_t& base, uint64_t delta) {
        if (delta == 1) {
            return base;
        } else if (delta % 2 == 0) {
            return as_handle(as_integer(base) + delta/2);
        } else {
            return as_handle(as_integer(base) - (delta-1)/2);
        }
    };

    auto edge_to_delta = [&](const edge_t& edge) {
        int64_t delta = as_integer(edge.second) - as_integer(edge.first);
        return (delta == 0 ? 1 : (delta > 0 ? 2*abs(delta) : 2*abs(delta)+1));
    };

    auto mask_edge = [&](const edge_t& edge) {
        uint64_t idx = masked_edges_bv.select1(as_integer(edge.first));
        masked_edges_iv.insert(idx+1, edge_to_delta(edge));
        masked_edges_bv.insert(idx+1, 0);
    };

    auto is_masked_edge = [&](const edge_t& edge) {
        uint64_t from_idx = as_integer(edge.first);
        uint64_t to_idx = as_integer(edge.second);
        uint64_t begin_edge = masked_edges_bv.select1(from_idx);
        uint64_t end_edge = masked_edges_bv.select1(from_idx+1);
        for (uint64_t i = begin_edge; i < end_edge; ++i) {
            if (as_integer(delta_to_handle(edge.first, masked_edges_iv.at(i))) == to_idx) {
                return true;
            }
        }
        return false;
    };
    
    // This (s) is our set of oriented nodes.
    // using a map instead of a set ensures a stable sort across different systems
    //map<handlegraph::nid_t, handle_t> s;
    dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,256,16> > s;
    uint64_t max_handle_rank = 0;
    g->for_each_handle([&](const handle_t& found) {
            max_handle_rank = std::max(max_handle_rank,
                                       number_bool_packing::unpack_number(found));
        });
    for (uint64_t i = 0; i <= max_handle_rank+1; ++i) {
        s.push_back(0);
        masked_edges_bv.push_back(1);
        masked_edges_bv.push_back(1);
        masked_edges_iv.push_back(0);
        masked_edges_iv.push_back(0);
    }

    // We find the head and tails, if there are any
    //std::vector<handle_t> heads{head_nodes(g)};
    // No need to fetch the tails since we don't use them
    
    // Maps from node ID to first orientation we suggested for it.
    dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,256,16> > seeds;
    dyn::hacked_vector seeds_rev;

    // Dump all the heads into the oriented set, rather than having them as
    // seeds. We will only go for cycle-breaking seeds when we run out of
    // heads. This is bad for contiguity/ordering consistency in cyclic
    // graphs and reversing graphs, but makes sure we work out to just
    // topological sort on DAGs. It mimics the effect we used to get when we
    // joined all the head nodes to a new root head node and seeded that. We
    // ignore tails since we only orient right from nodes we pick.
    if (use_heads) {
        for(const handle_t& head : head_nodes(g)) {
            s.set(number_bool_packing::unpack_number(head), 1);
        }
    } else if (use_tails) {
        for(const handle_t& tail : tail_nodes(g)) {
            s.set(number_bool_packing::unpack_number(tail), 1);
        }
    }

    // We will use an ordered map handles by ID for nodes we have not visited
    // yet. This ensures a consistent sort order across systems.
    dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,256,16> > unvisited;
    //map<handlegraph::nid_t, handle_t> unvisited;
    for (uint64_t i = 0; i <= max_handle_rank; ++i) {
        unvisited.push_back(0);
        seeds.push_back(0);
        seeds_rev.push_back(0);
    }
    g->for_each_handle([&](const handle_t& found) {
            uint64_t rank = number_bool_packing::unpack_number(found);
            unvisited.set(rank, !s.at(rank));
    });

    std::unique_ptr<progress_meter::ProgressMeter> progress;
    if (progress_reporting) {
        std::string banner = "[odgi::topological_order] sorting nodes:";
        progress = std::make_unique<progress_meter::ProgressMeter>(g->get_node_count(), banner);
    }

    while(unvisited.rank1(unvisited.size())!=0 || s.rank1(s.size())!=0) {

        // Put something in s. First go through seeds until we can find one
        // that's not already oriented.
        while(s.rank1(s.size())==0 && seeds.rank1(seeds.size())!=0) {
            // Look at the first seed
            //auto first_seed = (*seeds.begin()).second;
            uint64_t seed_rank = seeds.select1(0);
            handle_t first_seed = number_bool_packing::pack(seed_rank, seeds_rev.at(seed_rank));

            if(unvisited.at(number_bool_packing::unpack_number(first_seed))) {
                // We have an unvisited seed. Use it
#ifdef debug
#pragma omp critical (cerr)
                cerr << "Starting from seed " << g->get_id(first_seed) << " orientation " << g->get_is_reverse(first_seed) << endl;
#endif

                s.set(number_bool_packing::unpack_number(first_seed), 1);
                unvisited.set(number_bool_packing::unpack_number(first_seed), 0);
            }
            // Whether we used the seed or not, don't keep it around
            //seeds.erase(seeds.begin());
            seeds[seed_rank] = 0;
        }

        if(s.rank1(s.size())==0) {
            // If we couldn't find a seed, just grab any old node.
            // Since map order is stable across systems, we can take the first node by id and put it locally forward.
            /*
#ifdef debug
#pragma omp critical (cerr)
            cerr << "Starting from arbitrary node " << unvisited.select1(0) << " locally forward" << endl;
#endif
            */

            uint64_t h = unvisited.select1(0);
            s.set(h, 1);
            unvisited.set(h, 0);
        }

        while (s.rank1(s.size())!=0) {
            // Grab an oriented node
            uint64_t i = s.select1(0);
            s.set(i, 0);
            handle_t n = number_bool_packing::pack(i, false);
            // Emit it
            sorted.push_back(n);
            if (progress_reporting) {
                progress->increment(1);
            }
#ifdef debug
#pragma omp critical (cerr)
            cerr << "Using oriented node " << g->get_id(n) << " orientation " << g->get_is_reverse(n) << endl;
#endif

            // See if it has an edge from its start to the start of some node
            // where both were picked as places to break into cycles. A
            // reversing self loop on a cycle entry point is a special case of
            // this.
            g->follow_edges(n, true, [&](const handle_t& prev_node) {
                    if(!unvisited.at(number_bool_packing::unpack_number(prev_node))) {
                        // Look at the edge
                        auto edge = g->edge_handle(prev_node, n);
                        if (is_masked_edge(edge)) {
                            // We removed this edge, so skip it.
                            return;
                        }
                    
#ifdef debug
#pragma omp critical (cerr)
                        cerr << "\tHas left-side edge to cycle entry point " << g->get_id(prev_node)
                             << " orientation " << g->get_is_reverse(prev_node) << endl;
#endif

                        // Mask the edge
                        mask_edge(edge);
                    
#ifdef debug
#pragma omp critical (cerr)
                        cerr << "\t\tEdge: " << g->get_id(edge.first) << " " << g->get_is_reverse(edge.first)
                             << " -> " << g->get_id(edge.second) << " " << g->get_is_reverse(edge.second) << endl;
#endif
                    }
                });

            // All other connections and self loops are handled by looking off the right side.

            // See what all comes next, minus deleted edges.
            g->follow_edges(n, false, [&](const handle_t& next_node) {

                    uint64_t next_node_rank = number_bool_packing::unpack_number(next_node);

                    // Look at the edge
                    auto edge = g->edge_handle(n, next_node);
                    if (is_masked_edge(edge)) {
                        // We removed this edge, so skip it.
                        return;
                    }

#ifdef debug
#pragma omp critical (cerr)
                    cerr << next_node_rank << endl;
                    cerr << "\tHas edge to " << g->get_id(next_node) << " orientation " << g->get_is_reverse(next_node) << endl;
#endif

                    // Mask the edge connecting these nodes in this order and
                    // relative orientation, so we can't traverse it again

#ifdef debug
#pragma omp critical (cerr)
                    cerr << "\t\tEdge: " << g->get_id(edge.first) << " " << g->get_is_reverse(edge.first)
                         << " -> " << g->get_id(edge.second) << " " << g->get_is_reverse(edge.second) << endl;
#endif

                    // Mask the edge
                    mask_edge(edge);

                    if(unvisited.at(next_node_rank)) {
                        // We haven't already started here as an arbitrary cycle entry point

#ifdef debug
#pragma omp critical (cerr)
                        cerr << "\t\tAnd node hasn't been visited yet" << endl;
#endif

                        bool unmasked_incoming_edge = false;
                        g->follow_edges(next_node, true, [&](const handle_t& prev_node) {
                                // Get a handle for each incoming edge
                                auto prev_edge = g->edge_handle(prev_node, next_node);

                                if (!is_masked_edge(prev_edge)) {
                                    // We found such an edghe and can stop looking
                                    unmasked_incoming_edge = true;
                                    return false;
                                }
                                // Otherwise check all the edges on the left of this handle
                                return true;
                            });

                        if(!unmasked_incoming_edge) {

#ifdef debug
#pragma omp critical (cerr)
                            cerr << "\t\t\tIs last incoming edge" << endl;
#endif
                            // Keep this orientation and put it here
                            s.set(next_node_rank, 1);
                            // Remember that we've visited and oriented this node, so we
                            // don't need to use it as a seed.
                            unvisited.set(next_node_rank, 0);

                        } else if (!seeds[next_node_rank]) {
                            // We came to this node in this orientation; when we need a
                            // new node and orientation to start from (i.e. an entry
                            // point to the node's cycle), we might as well pick this
                            // one.
                            // Only take it if we don't already know of an orientation for this node.
                            seeds[next_node_rank] = 1;
                            seeds_rev[next_node_rank] = number_bool_packing::unpack_bit(next_node);

#ifdef debug
#pragma omp critical (cerr)
                            cerr << "\t\t\tSuggests seed " << g->get_id(next_node) << " orientation " << g->get_is_reverse(next_node) << endl;
#endif
                        }
                    } else {
#ifdef debug
#pragma omp critical (cerr)
                        cerr << "\t\tAnd node was already visited (to break a cycle)" << endl;
#endif
                    }
                });
        }
    }

    if (progress_reporting) {
        progress->finish();
    }

    assert(sorted.size() == g->get_node_count());

    // Send away our sorted ordering.
    return sorted;
}

std::vector<handle_t> two_way_topological_order(const HandleGraph* g) {
    // take the average assigned order for each handle
    hash_map<handle_t, uint64_t> avg_order;
    uint64_t i = 0;
    for (auto& handle : topological_order(g, true)) {
        avg_order[handle] = ++i;
    }
    i = 0;
    for (auto& handle : topological_order(g, false)) {
        avg_order[handle] = max(avg_order[handle], (++i));
    }
    std::vector<std::pair<handle_t, uint64_t>> order;
    for (auto& p : avg_order) {
        order.push_back(p);
    }
    std::sort(order.begin(), order.end(), [](const std::pair<handle_t, uint64_t>& a,
                                             const std::pair<handle_t, uint64_t>& b) {
                  return a.second < b.second;
              });
    std::vector<handle_t> result;
    for (auto& p : order) {
        result.push_back(p.first);
    }
    return result;
}

std::vector<handle_t> lazy_topological_order_internal(const HandleGraph* g, bool lazier) {
    
    // map that will contain the orientation and the in degree for each node
    hash_map<handle_t, int64_t> inward_degree;
    inward_degree.reserve(g->get_node_count());
    
    // stack for the traversal
    std::vector<handle_t> stack;
    
    if (lazier) {
        // take the locally forward orientation as a single stranded orientation
        g->for_each_handle([&](const handle_t& handle) {
            int64_t& degree = inward_degree[handle];
            g->follow_edges(handle, true, [&](const handle_t& ignored) {
                degree++;
            });
            // initialize the stack with head nodes
            if (degree == 0) {
                stack.emplace_back(handle);
            }
        });
    }
    else {
        // get an orientation over which we can consider the graph single stranded
        std::vector<handle_t> orientation = single_stranded_orientation(g);
        
        if (orientation.size() != g->get_node_count()) {
            cerr << "error:[algorithms] attempting to use lazy topological sort on unorientable graph" << endl;
            assert(false);
        }
        
        // compute the degrees by following the edges backward
        for (auto& handle : orientation) {
            int64_t& degree = inward_degree[handle];
            g->follow_edges(handle, true, [&](const handle_t& ignored) {
                degree++;
            });
            // initialize the stack with head nodes
            if (degree == 0) {
                stack.emplace_back(handle);
            }
        }
    }
    
    // the return value
    std::vector<handle_t> order;
    order.reserve(g->get_node_count());
    
    while (!stack.empty()) {
        // get a head node off the queue
        handle_t here = stack.back();
        stack.pop_back();
        
        // add it to the topological order
        order.push_back(here);
        
        // remove its outgoing edges
        g->follow_edges(here, false, [&](const handle_t& next) {
            
            auto iter = inward_degree.find(next);
            // we should never be able to reach the opposite orientation of a node
            assert(iter != inward_degree.end());
            // implicitly remove the edge
            iter->second--;
            if (iter->second == 0) {
                // after removing this edge, the node is now a head, add it to the queue
                stack.push_back(next);
            }
        });
    }
    
    if (order.size() != g->get_node_count()) {
        cerr << "error:[algorithms] lazy topological sort is invalid on non-DAG graph, cannot complete algorithm" << endl;
        assert(false);
    }
    
    return order;
}
    
    
std::vector<handle_t> lazy_topological_order(const HandleGraph* g) {
    return lazy_topological_order_internal(g, false);
}
    
std::vector<handle_t> lazier_topological_order(const HandleGraph* g) {
    return lazy_topological_order_internal(g, true);
}


void topological_sort(MutableHandleGraph& g, bool compact_ids) {
    g.apply_ordering(topological_order(&g), compact_ids);
}

std::vector<handle_t> breadth_first_topological_order(const HandleGraph& g, const uint64_t& chunk_size,
                                                      bool use_heads, bool use_tails) {

    // This (s) is our set of oriented nodes.
    //dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,256,16> > s;
    uint64_t min_handle_rank = 0;
    uint64_t max_handle_rank = 0;
    g.for_each_handle([&](const handle_t& found) {
                          uint64_t handle_rank = number_bool_packing::unpack_number(found);
                          min_handle_rank = std::min(min_handle_rank, handle_rank);
                          max_handle_rank = std::max(max_handle_rank, handle_rank);
                      });

    // Start with the heads of the graph.
    // We could also just use the first node of the graph.
    std::vector<handle_t> seeds;
    if (use_heads) {
        seeds = head_nodes(&g);
    } else if (use_tails) {
        seeds = tail_nodes(&g);
    } else {
        handle_t min_handle = number_bool_packing::pack(min_handle_rank, false);
        seeds = { min_handle };
    }

    // We need to keep track of the nodes we haven't visited to seed subsequent
    // runs of the BFS
    dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,256,16> > unvisited;
    for (uint64_t i = 0; i <= max_handle_rank; ++i) {
        unvisited.push_back(1);
    }
    /*
    g.for_each_handle([&](const handle_t& found) {
                          uint64_t rank = number_bool_packing::unpack_number(found);
                          unvisited.set(rank, 1);
                      });
    */

    uint64_t prev_max_root = 0;
    uint64_t prev_max_length = 0;
    
    std::vector<bfs_state_t> order_raw;
    while (unvisited.rank1(unvisited.size())!=0) {
        /*
        std::cerr << "unvisited size " << unvisited.rank1(unvisited.size()) << std::endl;
        for (uint64_t i = 0; i < unvisited.size(); ++i) {
            std::cerr << unvisited.at(i);
        }
        std::cerr << std::endl;
        */
        uint64_t seen_bp = 0;
        uint64_t curr_max_root = 0;
        uint64_t curr_max_length = 0;
        bfs(g,
            [&g,&order_raw,&unvisited,&seen_bp,
             &prev_max_root,&curr_max_root,
             &prev_max_length,&curr_max_length]
            (const handle_t& h, const uint64_t& r, const uint64_t& l, const uint64_t& d) {
                uint64_t i = number_bool_packing::unpack_number(h);
                order_raw.push_back({h, r+prev_max_root, l+prev_max_length});
                curr_max_root = std::max(r+prev_max_root, curr_max_root);
                curr_max_length = std::max(l+prev_max_length, curr_max_length);
                seen_bp += g.get_length(h);
                unvisited.set(i, 0);
            },
            [&unvisited](const handle_t& h) {
                uint64_t i = number_bool_packing::unpack_number(h);
                return unvisited.at(i)==0;
            },
            [](const handle_t& l, const handle_t& h) { return false; },
            [&seen_bp,&chunk_size]() { return seen_bp > chunk_size; },
            seeds,
            { },
            false); // don't use bidirectional search
        // get another seed
        prev_max_root = curr_max_root;
        prev_max_length = curr_max_length;
        if (unvisited.rank1(unvisited.size())!=0) {
            uint64_t i = unvisited.select1(0);
            handle_t h = number_bool_packing::pack(i, false);
            seeds = { h };
        }
    }
    //std::cerr << "order size " << order.size() << " graph size " << g.get_node_count() << std::endl;
    //assert(order.size() == g.get_node_count());
    std::sort(order_raw.begin(), order_raw.end(),
              [](const bfs_state_t& a,
                 const bfs_state_t& b) {
                  return a.root < b.root || a.root == b.root && a.length < b.length;
              });

    std::vector<handle_t> order;
    for (auto& o : order_raw) order.push_back(o.handle);

    return order;
}


std::vector<handle_t> depth_first_topological_order(const HandleGraph& g, const uint64_t& chunk_size,
                                                    bool use_heads, bool use_tails) {

    uint64_t min_handle_rank = 0;
    uint64_t max_handle_rank = 0;
    g.for_each_handle([&](const handle_t& found) {
                          uint64_t handle_rank = number_bool_packing::unpack_number(found);
                          min_handle_rank = std::min(min_handle_rank, handle_rank);
                          max_handle_rank = std::max(max_handle_rank, handle_rank);
                      });

    // Start with the heads of the graph.
    // We could also just use the first node of the graph.
    std::vector<handle_t> seeds;
    if (use_heads) {
        seeds = head_nodes(&g);
    } else if (use_tails) {
        seeds = tail_nodes(&g);
    } else {
        handle_t min_handle = number_bool_packing::pack(min_handle_rank, false);
        seeds = { min_handle };
    }

    // We need to keep track of the nodes we haven't visited to seed subsequent
    // runs of the BFS
    dyn::succinct_bitvector<dyn::spsi<dyn::packed_vector,256,16> > unvisited;
    for (uint64_t i = 0; i <= max_handle_rank; ++i) {
        unvisited.push_back(1);
    }
    /*
    g.for_each_handle([&](const handle_t& found) {
                          uint64_t rank = number_bool_packing::unpack_number(found);
                          unvisited.set(rank, 1);
                      });
    */
    std::vector<handle_t> order;
    while (unvisited.rank1(unvisited.size())!=0) {
        /*
        std::cerr << "unvisited size " << unvisited.rank1(unvisited.size()) << std::endl;
        for (uint64_t i = 0; i < unvisited.size(); ++i) {
            std::cerr << unvisited.at(i);
        }
        std::cerr << std::endl;
        */
        uint64_t bp_count = 0;
        dfs(g,
            [&g,&order,&unvisited,&bp_count](const handle_t& h) {
                uint64_t i = number_bool_packing::unpack_number(h);
                bp_count += g.get_length(h);
                order.push_back(h);
                unvisited.set(i, 0);
            },
            [](const handle_t& h) { },
            [&unvisited](const handle_t& h) {
                uint64_t i = number_bool_packing::unpack_number(h);
                return unvisited.at(i)==0;
            },
            [&bp_count,&chunk_size]() {
                //std::cerr << "bp_count " << bp_count << std::endl;
                return bp_count > chunk_size;
            },
            seeds);
        // get another seed
        if (unvisited.rank1(unvisited.size())!=0) {
            uint64_t i = unvisited.select1(0);
            handle_t h = number_bool_packing::pack(i, false);
            seeds = { h };
        }
    }
    //std::cerr << "order size " << order.size() << " graph size " << g.get_node_count() << std::endl;
    //assert(order.size() == g.get_node_count());

    return order;

}

}

}
