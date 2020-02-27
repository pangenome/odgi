#include "xp.hpp"

#include <arpa/inet.h>

namespace xp {

    ////////////////////////////////////////////////////////////////////////////
    // Here is XP
    ////////////////////////////////////////////////////////////////////////////

    XP::~XP(void) {
        // Clean up any created XGPaths
        while (!paths.empty()) {
            delete paths.back();
            paths.pop_back();
        }
    }

    /// build the graph from a graph handle
    void XP::from_handle_graph(const HandleGraph& graph) {
        // set up our enumerators
        auto for_each_sequence = [&](const std::function<void(const std::string& seq, const nid_t& node_id)>& lambda) {
            graph.for_each_handle([&](const handle_t& handle) {
                lambda(graph.get_sequence(handle), graph.get_id(handle));
            });
        };
        auto for_each_edge = [&](const std::function<void(const nid_t& from_id, const bool& from_rev,
                                                          const nid_t& to_id, const bool& to_rev)>& lambda) {
            graph.for_each_edge([&](const edge_t& edge) {
                lambda(graph.get_id(edge.first), graph.get_is_reverse(edge.first),
                       graph.get_id(edge.second), graph.get_is_reverse(edge.second));
            });
        };
        auto for_each_path_element = [&](const std::function<void(const std::string& path_name,
                                                                  const nid_t& node_id, const bool& is_rev,
                                                                  const std::string& cigar, const bool& is_empty, const bool& is_circular)>& lambda) {
            // no-op
        };
        from_enumerators(for_each_sequence, for_each_edge, for_each_path_element, false);
    }

    void XG::from_enumerators(const std::function<void(const std::function<void(const std::string& seq, const nid_t& node_id)>&)>& for_each_sequence,
                              const std::function<void(const std::function<void(const nid_t& from, const bool& from_rev,
                                                                                const nid_t& to, const bool& to_rev)>&)>& for_each_edge,
                              const std::function<void(const std::function<void(const std::string& path_name,
                                                                                const nid_t& node_id, const bool& is_rev,
                                                                                const std::string& cigar, const bool& is_empty,
                                                                                const bool& is_circular)>&)>& for_each_path_element,
                              bool validate, std::string basename) {

        if (basename.empty()) {
            basename = temp_file::create();
        }
        node_count = 0;
        seq_length = 0;
        edge_count = 0;
        path_count = 0;
        min_id = std::numeric_limits<int64_t>::max();
        max_id = 0;
        // get information about graph size and id ranges
#ifdef VERBOSE_DEBUG
        cerr << "computing graph sequence length and node count" << endl;
#endif
        for_each_sequence([&](const std::string& seq, const nid_t& id) {
            // min id starts at 0
            min_id = std::min(min_id, id);
            max_id = std::max(max_id, id);
            seq_length += seq.size();
            ++node_count;
        });
#ifdef VERBOSE_DEBUG
        cerr << "counting edges" << endl;
#endif
        // edge count
        for_each_edge([&](const nid_t& from_id, const bool& from_rev, const nid_t& to_id, const bool& to_rev) {
            ++edge_count;
        });
        // path count
        std::string pname;
#ifdef VERBOSE_DEBUG
        cerr << "counting paths" << endl;
#endif
        for_each_path_element([&](const std::string& path_name, const nid_t& node_id, const bool& is_rev, const std::string& cigar, const bool& is_empty, const bool& is_circular) {
            if (path_name != pname) {
                ++path_count;
            }
            pname = path_name;
        });
#ifdef VERBOSE_DEBUG
        std::cerr << "graph has " << seq_length << "bp in sequence, "
              << node_count << " nodes, "
              << edge_count << " edges, and "
              << path_count << " paths" << std::endl
              << "node ids run from " << min_id << " to " << max_id << std::endl;
#endif

        // FIXME CUT OFF FROM HERE
        /**
        // set up our compressed representation
        sdsl::int_vector<> i_iv;
        sdsl::util::assign(s_iv, sdsl::int_vector<>(seq_length, 0, 3));
        sdsl::util::assign(s_bv, sdsl::bit_vector(seq_length+1));
        sdsl::util::assign(i_iv, sdsl::int_vector<>(node_count));
        sdsl::util::assign(r_iv, sdsl::int_vector<>(max_id-min_id+1)); // note: possibly discontinuous

        // for each node in the sequence
        // concatenate the labels into the s_iv
#ifdef VERBOSE_DEBUG
        cerr << "storing node labels" << endl;
#endif
        size_t r = 1;
        // first make i_iv and r_iv
        for_each_sequence([&](const std::string& seq, const nid_t& id) {
            i_iv[r-1] = id;
            // store ids to rank mapping
            r_iv[id-min_id] = r;
            ++r;
        });
        sdsl::util::bit_compress(i_iv);
        sdsl::util::bit_compress(r_iv);

        // then make s_bv and s_iv
        size_t j = 0;
        for_each_sequence([&](const std::string& seq, const nid_t& id) {
            //size_t i = r_iv[id-min_id]-1;
            s_bv[j] = 1; // record node start
            for (auto c : seq) {
                s_iv[j++] = dna3bit(c); // store sequence
            }
        });
        s_bv[seq_length] = 1;

        // to label the paths we'll need to compress and index our vectors
        sdsl::util::bit_compress(s_iv);
        sdsl::util::assign(s_bv_rank, sdsl::rank_support_v<1>(&s_bv));
        sdsl::util::assign(s_bv_select, sdsl::bit_vector::select_1_type(&s_bv));

        // now that we've set up our sequence indexes, we can build the locally traversable graph storage

        auto temp_get_handle = [&](const nid_t& id, bool orientation) {
            uint64_t handle_rank = r_iv[id-min_id];
            return number_bool_packing::pack(handle_rank, orientation);
        };
        auto temp_node_size = [&](const nid_t& id) {
            uint64_t handle_rank = r_iv[id-min_id];
            return s_bv_select(handle_rank+1)-s_bv_select(handle_rank);
        };
        auto temp_get_id = [&](const handle_t& h) {
            return i_iv[number_bool_packing::unpack_number(h)-1];
        };

#ifdef VERBOSE_DEBUG
        cerr << "collecting edges " << endl;
#endif

        // first, we need to collect the edges for each node
        // we use the mmmultimap here to reduce in-memory costs to a minimum
        std::string edge_f_t_idx = basename + ".from_to.mm";
        std::string edge_t_f_idx = basename + ".to_from.mm";
        auto edge_from_to_mm = std::make_unique<mmmulti::map<uint64_t, uint64_t>>(edge_f_t_idx);
        auto edge_to_from_mm = std::make_unique<mmmulti::map<uint64_t, uint64_t>>(edge_t_f_idx);
        for_each_edge([&](const nid_t& from_id, const bool& from_rev, const nid_t& to_id, const bool& to_rev) {
            handle_t from_handle = temp_get_handle(from_id, from_rev);
            handle_t to_handle = temp_get_handle(to_id, to_rev);
            edge_from_to_mm->append(as_integer(from_handle), as_integer(to_handle));
            edge_to_from_mm->append(as_integer(to_handle), as_integer(from_handle));
        });
        handle_t max_handle = number_bool_packing::pack(r_iv.size(), true);
        edge_from_to_mm->index(as_integer(max_handle));
        edge_to_from_mm->index(as_integer(max_handle));

        // calculate g_iv size
        size_t g_iv_size =
                node_count * G_NODE_HEADER_LENGTH // record headers
                + edge_count * 2 * G_EDGE_LENGTH; // edges (stored twice)
        sdsl::util::assign(g_iv, sdsl::int_vector<>(g_iv_size));
        sdsl::util::assign(g_bv, sdsl::bit_vector(g_iv_size));

#ifdef VERBOSE_DEBUG
        cerr << "computing graph vector " << endl;
#endif

        int64_t g = 0; // pointer into g_iv and g_bv
        for (int64_t i = 0; i < node_count; ++i) {
            nid_t id = i_iv[i];
#ifdef VERBOSE_DEBUG
            if (i % 1000 == 0) cerr << i << " of " << node_count << " ~ " << (float)i/(float)node_count * 100 << "%" << "\r";
#endif
            handle_t handle = temp_get_handle(id, false);
            //std::cerr << "id " << id << std::endl;
            g_bv[g] = 1; // mark record start for later query
            g_iv[g++] = id;
            g_iv[g++] = node_vector_offset(id);
            g_iv[g++] = temp_node_size(id);
            size_t to_edge_count = 0;
            size_t from_edge_count = 0;
            size_t to_edge_count_idx = g++;
            size_t from_edge_count_idx = g++;
            // write the edges in id-based format
            // we will next convert these into relative format
            for (auto orientation : { false, true }) {
                handle_t to = temp_get_handle(id, orientation);
                //std::cerr << "looking at to handle " << number_bool_packing::unpack_number(to) << ":" << number_bool_packing::unpack_bit(to) << std::endl;
                edge_to_from_mm->for_unique_values_of(as_integer(to), [&](const uint64_t& _from) {
                    handle_t from = as_handle(_from);
                    //std::cerr << "edge to " << number_bool_packing::unpack_number(from) << ":" << number_bool_packing::unpack_bit(from)
                    //<< " -> " << number_bool_packing::unpack_number(to) << ":" << number_bool_packing::unpack_bit(to) << std::endl;
                    g_iv[g++] = temp_get_id(from);
                    g_iv[g++] = edge_type(from, to);
                    ++to_edge_count;
                });
            }
            g_iv[to_edge_count_idx] = to_edge_count;
            for (auto orientation : { false, true }) {
                handle_t from = temp_get_handle(id, orientation);
                //std::cerr << "looking at from handle " << number_bool_packing::unpack_number(from) << ":" << number_bool_packing::unpack_bit(from) << std::endl;
                edge_from_to_mm->for_unique_values_of(as_integer(from), [&](const uint64_t& _to) {
                    handle_t to = as_handle(_to);
                    //std::cerr << "edge from " << number_bool_packing::unpack_number(from) << ":" << number_bool_packing::unpack_bit(from)
                    //<< " -> " << number_bool_packing::unpack_number(to) << ":" << number_bool_packing::unpack_bit(to) << std::endl;
                    g_iv[g++] = temp_get_id(to);
                    g_iv[g++] = edge_type(from, to);
                    ++from_edge_count;
                });
            }
            g_iv[from_edge_count_idx] = from_edge_count;
        }

#ifdef VERBOSE_DEBUG
        std::cerr << node_count << " of " << node_count << " ~ 100.0000%" << std::endl;
#endif

        // cleanup our mmmultimap
        edge_from_to_mm.reset();
        edge_to_from_mm.reset();
        std::remove(edge_f_t_idx.c_str());
        std::remove(edge_t_f_idx.c_str());

        // set up rank and select supports on g_bv so we can locate nodes in g_iv
        sdsl::util::assign(g_bv_rank, sdsl::rank_support_v<1>(&g_bv));
        sdsl::util::assign(g_bv_select, sdsl::bit_vector::select_1_type(&g_bv));

#ifdef VERBOSE_DEBUG
        cerr << "making graph vector relativistic " << endl;
#endif

        // convert the edges in g_iv to relativistic form
        for (int64_t i = 0; i < node_count; ++i) {
            int64_t id = i_iv[i];
            // find the start of the node's record in g_iv
            int64_t g = g_bv_select(id_to_rank(id));
            // get to the edges to
            int edges_to_count = g_iv[g+G_NODE_TO_COUNT_OFFSET];
            int edges_from_count = g_iv[g+G_NODE_FROM_COUNT_OFFSET];
            int64_t t = g + G_NODE_HEADER_LENGTH;
            int64_t f = g + G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * edges_to_count;
            for (int64_t j = t; j < f; ) {
                g_iv[j] = g_bv_select(id_to_rank(g_iv[j])) - g;
                j += 2;
            }
            for (int64_t j = f; j < f + G_EDGE_LENGTH * edges_from_count; ) {
                g_iv[j] = g_bv_select(id_to_rank(g_iv[j])) - g;
                j += 2;
            }
        }
        sdsl::util::clear(i_iv);
        sdsl::util::bit_compress(g_iv);
        **/
        // FIXME END OF CUT OFF

#ifdef VERBOSE_DEBUG
        cerr << "storing paths" << endl;
#endif
        // paths
        std::string path_names;

        std::string curr_path_name;
        std::vector<handle_t> curr_path_steps;
        size_t curr_node_count = 0;
        bool curr_is_circular = false; // TODO, use TP:Z:circular tag... we'll have to fish this out of the file
        uint64_t p = 0;

        auto build_accumulated_path = [&](void) {
            // only build if we had a path to build
#ifdef VERBOSE_DEBUG
            if (++p % 100 == 0) std::cerr << p << " of " << path_count << " ~ " << (float)p/(float)path_count * 100 << "%" << "\r";
#endif

#ifdef debug_path_index
            std::cerr << "Creating XPPath for path " << curr_path_name << " with visits:" << std::endl;
        for (handle_t visiting : curr_path_steps) {
             std::cerr << "\tHandle for " << get_id(visiting) << (get_is_reverse(visiting) ? "-" : "+") << " with value "
                    << as_integer(visiting) << std::endl;
        }
#endif
            size_t unique_member_count = 0;
            path_names += start_marker + curr_path_name + end_marker;
            XPPath* path = new XPPath(curr_path_name, curr_path_steps,
                                      curr_is_circular,
                                      *this);
            paths.push_back(path);

#ifdef debug_path_index
            std::cerr << "paths[" << paths.size() - 1 << "] = " << curr_path_name << " @ " << path << std::endl;
#endif
        };

        // todo ... is it circular?
        // might make sense to scan the file for this
        bool has_path = false;
        for_each_path_element([&](const std::string& path_name, const nid_t& node_id, const bool& is_rev, const std::string& cigar, const bool& is_empty, const bool& is_circular) {
            if (path_name != curr_path_name && !curr_path_name.empty()) {
                // build the last path we've accumulated
                build_accumulated_path();
                curr_path_steps.clear();
                curr_is_circular = false;
            }
            curr_path_name = path_name;
            if (!is_empty) {
                handle_t visiting = get_handle(node_id, is_rev);
#ifdef debug_path_index
                std::cerr << "Adding handle for " << node_id << (is_rev ? "-" : "+") << " with value "
                    << as_integer(visiting) << " to path " << path_name << std::endl;
                std::cerr << "Node length is " << get_length(visiting) << " bp" << std::endl;
#endif
                curr_path_steps.push_back(visiting);
            }
            curr_is_circular = is_circular;
            has_path = true;
        });
        // build the last path
        if (has_path) {
            build_accumulated_path();
        }
        curr_path_steps.clear();
        curr_is_circular = false;
#ifdef VERBOSE_DEBUG
        std::cerr << path_count << " of " << path_count << " ~ 100.0000%" << std::endl;
#endif

        //std::cerr << "path names " << path_names << std::endl;

        // handle path names
        sdsl::util::assign(pn_iv, sdsl::int_vector<>(path_names.size()));
        sdsl::util::assign(pn_bv, sdsl::bit_vector(path_names.size()));
        // now record path name starts
        for (size_t i = 0; i < path_names.size(); ++i) {
            pn_iv[i] = path_names[i];
            if (path_names[i] == start_marker) {
                pn_bv[i] = 1; // register name start
            }
        }
        sdsl::util::assign(pn_bv_rank, sdsl::rank_support_v<1>(&pn_bv));
        sdsl::util::assign(pn_bv_select, sdsl::bit_vector::select_1_type(&pn_bv));

        // is this file removed by construct?
        string path_name_file = basename + ".pathnames.iv";
        sdsl::store_to_file((const char*)path_names.c_str(), path_name_file);
        sdsl::construct(pn_csa, path_name_file, 1);

        // FIXME CUT OFF FROM HERE
#ifdef VERBOSE_DEBUG
        cerr << "computing node to path membership" << endl;
#endif

        // create the node-to-path indexes
        index_node_to_path(basename);

        // validate the graph
        if (validate) {
            // do we get the correct sequences when looking up handles by id?
            for_each_sequence([&](const std::string& seq, const nid_t& id) {
                handle_t handle = get_handle(id);
                if (seq != get_sequence(handle)) {
                    std::cerr << "mismatch in handle sequence for " << id << std::endl;
                    exit(1);
                }
                if (get_id(handle) != id) {
                    std::cerr << "mismatch in id for " << id << std::endl;
                    exit(1);
                }
            });
            // do we have the correct set of edges?
            for_each_edge([&](const nid_t& from_id, const bool& from_rev, const nid_t& to_id, const bool& to_rev) {
                handle_t from_handle = get_handle(from_id, from_rev);
                handle_t to_handle = get_handle(to_id, to_rev);
                bool seen_to = false;
                follow_edges(from_handle, false, [&](const handle_t& h) {
                    //std::cerr << "fwd I see edge " << get_id(from_handle) << ":" << get_is_reverse(from_handle) << " -> " << get_id(h) << ":" << get_is_reverse(h) << std::endl;
                    //std::cerr << as_integer(h) << " ==? " << as_integer(to_handle) << std::endl;
                    if (h == to_handle) {
                        seen_to = true;
                    }
                });
                bool seen_from = false;
                follow_edges(to_handle, true, [&](const handle_t& h) {
                    //std::cerr << "rev I see edge " << get_id(h) << ":" << get_is_reverse(h) << " -> " << get_id(to_handle) << ":" << get_is_reverse(to_handle) << std::endl;
                    //std::cerr << as_integer(h) << " ==? " << as_integer(from_handle) << std::endl;
                    if (h == from_handle) {
                        seen_from = true;
                    }
                });
                if (!seen_to) {
                    std::cerr << "can't find to edge for " << get_id(from_handle) << ":" << get_is_reverse(from_handle)
                              << " -> " << get_id(to_handle) << ":" << get_is_reverse(to_handle) << std::endl;
                    exit(1);
                }
                if (!seen_from) {
                    std::cerr << "can't find from edge for " << get_id(from_handle) << ":" << get_is_reverse(from_handle)
                              << " -> " << get_id(to_handle) << ":" << get_is_reverse(to_handle) << std::endl;
                    exit(1);
                }
            });
            // do our stored paths match those in the input?

            std::string curr_path_name;
            std::vector<handle_t> curr_path_steps;
            size_t curr_node_count = 0;
            bool curr_is_circular = false; // TODO, use TP:Z:circular tag... we'll have to fish this out of the file
            uint64_t p_handle = 0;
            auto check_accumulated_path = [&](void) {
                ++p_handle;
                // only check if we had a path to build
                if (curr_path_steps.empty()) return;
                auto& path = *paths[p_handle-1];
                // check that the path name is correct
                if (get_path_name(as_path_handle(p_handle)) != curr_path_name) {
                    std::cerr << "path name mismatch " << get_path_name(as_path_handle(p_handle)) << " != " << curr_path_name << std::endl;
                    exit(1);
                }
                // check that the path handles are correct
                for (uint64_t i = 0; i < curr_path_steps.size(); ++i) {
                    if (path.handle(i) != curr_path_steps[i]) {
                        std::cerr << "handle mismatch " << get_id(path.handle(i)) << " != " << get_id(curr_path_steps[i]) << " in path " << curr_path_name << std::endl;
                        exit(1);
                    }
                }
                // check that the path positions are correct
                uint64_t pos = 0;
                path_handle_t path_handle = as_path_handle(p_handle);
                //std::cerr << "on path " << curr_path_name << std::endl;
                //const std::function<bool(const step_handle_t&, const bool&, const uint64_t&)>& iteratee) const;
                for (uint64_t i = 0; i < curr_path_steps.size(); ++i) {
                    handle_t handle = curr_path_steps[i];
                    //std::cerr << "looking at node " << get_id(handle) << " on path " << curr_path_name << std::endl;
                    uint64_t handle_length = get_length(handle);
                    for (uint64_t j = 0; j < handle_length; ++j) {
                        handle_t handle_at = path.handle_at_position(pos+j);
                        if (handle != handle_at) {
                            std::cerr << "handle at position mismatch " << get_id(handle) << " != " << get_id(handle_at)
                                      << " in path " << curr_path_name << " at position " << pos+j << std::endl;
                            exit(1);
                        }
                    }
                    bool path_seen = false;
                    auto check_pos_index = [&](const step_handle_t& step, const bool& is_rev, const uint64_t& pos) {
                        // check that the step handle is the same as this handle
                        //std::cerr << "checking " << get_id(handle) << " on " << curr_path_name << std::endl;
                        path_handle_t path_handle_of = get_path_handle_of_step(step);
                        if (path_handle == path_handle_of && as_integers(step)[1] == i) {
                            //std::cerr << "found matching step handle" << std::endl;
                            handle_t oriented_handle = !get_is_reverse(handle) && is_rev ? flip(handle) : handle;
                            handle_t handle_at = get_handle_of_step(step);
                            if (oriented_handle != handle_at) {
                                std::cerr << "handle at step mismatch " << get_id(oriented_handle) << " != " << get_id(handle_at)
                                          << " in path " << curr_path_name << std::endl;
                                std::cerr << "handle " << as_integer(handle_at) << " != " << as_integer(oriented_handle) << std::endl;
                                exit(1);
                            }
                            handle_t handle_at_pos = path.handle_at_position(pos);
                            if (handle_at != handle_at_pos) {
                                std::cerr << "handle at position mismatch " << get_id(handle_at) << " != " << get_id(handle_at_pos)
                                          << " in path " << curr_path_name << " at position " << pos << std::endl;
                                exit(1);
                            }
                            path_seen = true;
                            return false;
                        }
                        return true;
                    };
                    // verify that we have the right position in the node to path position index
                    for_each_step_position_on_handle(handle, check_pos_index);
                    if (!path_seen) {
                        std::cerr << "didn't find path " << curr_path_name << " in reverse index for " << get_id(handle) << std::endl;
                        exit(1);
                    }
                    pos += get_length(handle);
                }
            };
            for_each_path_element([&](const std::string& path_name, const nid_t& node_id, const bool& is_rev, const std::string& cigar, const bool& is_empty, const bool& is_circular) {
                if (path_name != curr_path_name && !curr_path_name.empty()) {
                    // check the last path we've accumulated
                    check_accumulated_path();
                    curr_path_steps.clear();
                }
                curr_path_name = path_name;
                if (!is_empty) {
                    curr_path_steps.push_back(get_handle(node_id, is_rev));
                }
            });
            // check the last path
            check_accumulated_path();
            curr_path_steps.clear();
        }

//#define DEBUG_CONSTRUCTION
#ifdef DEBUG_CONSTRUCTION
        cerr << "|g_iv| = " << size_in_mega_bytes(g_iv) << endl;
    cerr << "|g_bv| = " << size_in_mega_bytes(g_bv) << endl;
    cerr << "|s_iv| = " << size_in_mega_bytes(s_iv) << endl;

    cerr << "|r_iv| = " << size_in_mega_bytes(r_iv) << endl;

    cerr << "|s_bv| = " << size_in_mega_bytes(s_bv) << endl;

    long double paths_mb_size = 0;
    cerr << "|pn_iv| = " << size_in_mega_bytes(pn_iv) << endl;
    paths_mb_size += size_in_mega_bytes(pn_iv);
    cerr << "|pn_csa| = " << size_in_mega_bytes(pn_csa) << endl;
    paths_mb_size += size_in_mega_bytes(pn_csa);
    cerr << "|pn_bv| = " << size_in_mega_bytes(pn_bv) << endl;
    paths_mb_size += size_in_mega_bytes(pn_bv);
    paths_mb_size += size_in_mega_bytes(pn_bv_rank);
    paths_mb_size += size_in_mega_bytes(pn_bv_select);
    paths_mb_size += size_in_mega_bytes(pi_iv);
    cerr << "|np_iv| = " << size_in_mega_bytes(np_iv) << endl;
    paths_mb_size += size_in_mega_bytes(np_iv);
    cerr << "|np_bv| = " << size_in_mega_bytes(np_bv) << endl;
    paths_mb_size += size_in_mega_bytes(np_bv);
    paths_mb_size += size_in_mega_bytes(np_bv_select);
    cerr << "total paths size " << paths_mb_size << endl;

    float path_ids_mb_size=0;
    float path_dir_mb_size=0;
    float path_pos_mb_size=0;
    float path_ranks_mb_size=0;
    float path_offsets_mb_size=0;
    for (size_t i = 0; i < paths.size(); i++) {
        // Go through paths by number, so we can determine rank
        XGPath* path = paths[i];
        path_ids_mb_size += size_in_mega_bytes(path->ids);
        path_dir_mb_size += size_in_mega_bytes(path->directions);
        path_pos_mb_size += size_in_mega_bytes(path->positions);
        path_ranks_mb_size += size_in_mega_bytes(path->ranks);
        path_offsets_mb_size += size_in_mega_bytes(path->offsets);
    }
    cerr << "path ids size " << path_ids_mb_size << endl;
    cerr << "path directions size " << path_dir_mb_size << endl;
    cerr << "path positions size " << path_pos_mb_size << endl;
    cerr << "path ranks size " << path_ranks_mb_size << endl;
    cerr << "path offsets size " << path_offsets_mb_size << endl;

#endif

        bool print_graph = false;
        if (print_graph) {
            cerr << "printing graph" << endl;
            // we have to print the relativistic graph manually because the default sdsl printer assumes unsigned integers are stored in it
            for (size_t i = 0; i < g_iv.size(); ++i) {
                cerr << (int64_t)g_iv[i] << " ";
            } cerr << endl;
            for (int64_t i = 0; i < node_count; ++i) {
                int64_t id = rank_to_id(i+1);
                // find the start of the node's record in g_iv
                int64_t g = g_bv_select(id_to_rank(id));
                // get to the edges to
                int edges_to_count = g_iv[g+G_NODE_TO_COUNT_OFFSET];
                int edges_from_count = g_iv[g+G_NODE_FROM_COUNT_OFFSET];
                int sequence_size = g_iv[g+G_NODE_LENGTH_OFFSET];
                size_t seq_start = g_iv[g+G_NODE_SEQ_START_OFFSET];
                cerr << id << " ";
                for (int64_t j = seq_start; j < seq_start+sequence_size; ++j) {
                    cerr << revdna3bit(s_iv[j]);
                } cerr << " : ";
                int64_t t = g + G_NODE_HEADER_LENGTH;
                int64_t f = g + G_NODE_HEADER_LENGTH + G_EDGE_LENGTH * edges_to_count;
                cerr << " from ";
                for (int64_t j = t; j < f; ) {
                    cerr << rank_to_id(g_bv_rank(g+g_iv[j])+1) << " ";
                    j += 2;
                }
                cerr << " to ";
                for (int64_t j = f; j < f + G_EDGE_LENGTH * edges_from_count; ) {
                    cerr << rank_to_id(g_bv_rank(g+g_iv[j])+1) << " ";
                    j += 2;
                }
                cerr << endl;
            }
            cerr << s_iv << endl;
            for (size_t i = 0; i < s_iv.size(); ++i) {
                cerr << revdna3bit(s_iv[i]);
            } cerr << endl;
            cerr << s_bv << endl;
            cerr << "paths (" << paths.size() << ")" << endl;
            for (size_t i = 0; i < paths.size(); i++) {
                // Go through paths by number, so we can determine rank
                XGPath* path = paths[i];
                cerr << get_path_name(as_path_handle(i + 1)) << endl;
                // manually print IDs because simplified wavelet tree doesn't support ostream for some reason
                for (size_t j = 0; j + 1 < path->handles.size(); j++) {
                    cerr << get_id(path->handle(j)) << " ";
                }
                if (path->handles.size() > 0) {
                    cerr << get_id(path->handle(path->handles.size() - 1));
                }
                cerr << endl;
                cerr << path->offsets << endl;
            }
            cerr << np_bv << endl;
            cerr << np_iv << endl;
            cerr << nx_iv << endl;

        }

    }

    // FIXME This was directly copied from xg.cpp
    uint32_t XP::get_magic_number(void) const {
        return 4143290017ul;
    }

    void XG::serialize_members(ostream& out) const {
        serialize_and_measure(out);
    }

    // TODO Clean this up.
    size_t XG::serialize_and_measure(ostream& out, sdsl::structure_tree_node* s, std::string name) const {

        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
        size_t written = 0;

        // Do the magic number
        out << "XP";
        written += 2;

        ////////////////////////////////////////////////////////////////////////
        // DO NOT CHANGE THIS CODE without creating a new XG version:
        // 1. Increment CURRENT_VERSION to a new integer.
        // 2. Add your new serialization code.
        // 3. Add a new case for your new version to XG::load()
        ////////////////////////////////////////////////////////////////////////

        written += sdsl::write_member(s_iv.size(), out, child, "sequence_length");
        written += sdsl::write_member(node_count, out, child, "node_count");
        written += sdsl::write_member(edge_count, out, child, "edge_count");
        written += sdsl::write_member(path_count, out, child, "path_count");
        written += sdsl::write_member(min_id, out, child, "min_id");
        written += sdsl::write_member(max_id, out, child, "max_id");

        written += r_iv.serialize(out, child, "rank_id_vector");

        written += g_iv.serialize(out, child, "graph_vector");
        written += g_bv.serialize(out, child, "graph_bit_vector");
        written += g_bv_rank.serialize(out, child, "graph_bit_vector_rank");
        written += g_bv_select.serialize(out, child, "graph_bit_vector_select");

        written += s_iv.serialize(out, child, "seq_vector");
        written += s_bv.serialize(out, child, "seq_node_starts");
        written += s_bv_rank.serialize(out, child, "seq_node_starts_rank");
        written += s_bv_select.serialize(out, child, "seq_node_starts_select");

        // Treat the paths as their own node
        size_t paths_written = 0;
        auto paths_child = sdsl::structure_tree::add_child(child, "paths", sdsl::util::class_name(*this));

        paths_written += pn_iv.serialize(out, paths_child, "path_names");
        paths_written += pn_csa.serialize(out, paths_child, "path_names_csa");
        paths_written += pn_bv.serialize(out, paths_child, "path_names_starts");
        paths_written += pn_bv_rank.serialize(out, paths_child, "path_names_starts_rank");
        paths_written += pn_bv_select.serialize(out, paths_child, "path_names_starts_select");
        paths_written += pi_iv.serialize(out, paths_child, "path_ids");
        // TODO: Path count is written twice (once from paths.size() and once earlier from path_count)
        // We should remove one and cut a new xg version
        paths_written += sdsl::write_member(paths.size(), out, paths_child, "path_count");
        for (size_t i = 0; i < paths.size(); i++) {
            XGPath* path = paths[i];
            paths_written += path->serialize(out, paths_child, "path:" + get_path_name(as_path_handle(i+1)));
        }

        paths_written += np_bv.serialize(out, paths_child, "node_path_mapping_starts");
        paths_written += np_bv_select.serialize(out, paths_child, "node_path_mapping_starts_select");
        paths_written += np_iv.serialize(out, paths_child, "node_path_mapping");
        paths_written += nr_iv.serialize(out, paths_child, "node_path_rank");
        paths_written += nx_iv.serialize(out, paths_child, "node_path_position");

        sdsl::structure_tree::add_size(paths_child, paths_written);
        written += paths_written;

        sdsl::structure_tree::add_size(child, written);
        return written;

    }

    void XP::deserialize_members(std::istream& in) {
        // simple alias to match an external interface
        load(in);
    }

    // FIXME This does not work as it is.
    void XG::load(std::istream& in) {

        if (!in.good()) {
            throw XGFormatError("Index file does not exist or index stream cannot be read");
        }

        // Version 0 is the last XG format without an explicit version specifier.
        // If we find a version specifier we will up this.
        uint32_t file_version = 0;

        // We need to look for the magic value
        char buffer;
        in.get(buffer);
        if (buffer == 'X') {
            in.get(buffer);
            if (buffer == 'G') {
                // We found the magic value!

                // Don't put it back, but the next 4 bytes are a version number.
                in.read((char*) &file_version, sizeof(file_version));
                // Make sure to convert from network to host byte order
                file_version = ntohl(file_version);

            } else {
                // Put back both characters
                in.unget();
                in.unget();
            }
        } else {
            // Put back the one character
            in.unget();
        }

        if (file_version > CURRENT_VERSION) {
            // This XG file is from the distant future.
            throw XGFormatError("XG index file version " + std::to_string(file_version) +
                                " is too new to interpret (max version = " + std::to_string(CURRENT_VERSION) + ")");
        }

        try {

            ////////////////////////////////////////////////////////////////////////
            // DO NOT CHANGE THIS CODE without creating a new XG version:
            // 1. Increment OUTPUT_VERSION to a new integer.
            // 2. Change the serialization code.
            // 3. Add a new case here (or enhance an existing case) with new deserialization code.
            ////////////////////////////////////////////////////////////////////////
            switch (file_version) {

                case 5: // Fall through
                case 6:
                case 7:
                case 8:
                case 9:
                case 10:
                case 11:
                case 12:
                    std::cerr << "warning:[XG] Loading an out-of-date XG format."
                              << "For better performance over repeated loads, consider recreating this XG index." << std::endl;
                    // Fall through
                case 13:
                {
                    sdsl::read_member(seq_length, in);
                    sdsl::read_member(node_count, in);
                    sdsl::read_member(edge_count, in);
                    sdsl::read_member(path_count, in);
                    size_t entity_count = node_count + edge_count;
                    //cerr << sequence_length << ", " << node_count << ", " << edge_count << endl;
                    sdsl::read_member(min_id, in);
                    sdsl::read_member(max_id, in);

                    if (file_version <= 8) {
                        // Load the old id int vector to skip
                        sdsl::int_vector<> i_iv;
                        i_iv.load(in);
                    }
                    r_iv.load(in);

                    g_iv.load(in);
                    g_bv.load(in);
                    g_bv_rank.load(in, &g_bv);
                    g_bv_select.load(in, &g_bv);

                    s_iv.load(in);
                    s_bv.load(in);
                    s_bv_rank.load(in, &s_bv);
                    s_bv_select.load(in, &s_bv);

                    if (file_version <= 11) {
                        // Skip over gPBWT thread names
                        {
                            sdsl::csa_bitcompressed<> tn_csa;
                            tn_csa.load(in);
                        }
                        {
                            sdsl::sd_vector<> tn_cbv;
                            tn_cbv.load(in);
                            sdsl::sd_vector<>::rank_1_type tn_cbv_rank;
                            tn_cbv_rank.load(in, &tn_cbv);
                            sdsl::sd_vector<>::select_1_type tn_cbv_select;
                            tn_cbv_select.load(in, &tn_cbv);
                        }
                    }

                    if (file_version >= 7 && file_version <= 11) {
                        // There is a single haplotype count here for all components
                        // We ignore this part of the gPBWT
                        size_t haplotype_count;
                        sdsl::read_member(haplotype_count, in);
                    }

                    if (file_version <= 11) {
                        // Skip thread positions in gPBWT
                        {
                            sdsl::vlc_vector<> tin_civ;
                            tin_civ.load(in);
                        }
                        {
                            sdsl::vlc_vector<> tio_civ;
                            tio_civ.load(in);
                        }
                        {
                            sdsl::wt_int<> side_thread_wt;
                            side_thread_wt.load(in);
                        }
                    }

                    pn_iv.load(in);
                    pn_csa.load(in);
                    pn_bv.load(in);
                    pn_bv_rank.load(in, &pn_bv);
                    pn_bv_select.load(in, &pn_bv);
                    pi_iv.load(in);
                    sdsl::read_member(path_count, in);
                    for (size_t i = 0; i < path_count; ++i) {
                        auto path = new XGPath;
                        // Load the path, giving it the file version and a
                        // rank-to-ID comversion function for format upgrade
                        // purposes.
                        if (file_version <= 12) {
                            path->load_from_old_version(in, file_version, *this);
                        }
                        else {
                            path->load(in);
                        }
                        paths.push_back(path);
                    }

                    if (file_version <= 12) {
                        // skip over the old node-to-path indexes
                        {
                            sdsl::int_vector<> old_np_iv;
                            old_np_iv.load(in);
                        }
                        {
                            sdsl::bit_vector old_np_bv;
                            old_np_bv.load(in);
                            sdsl::rank_support_v<1> old_np_bv_rank;
                            old_np_bv_rank.load(in, &np_bv);
                            sdsl::bit_vector::select_1_type old_np_bv_select;
                            old_np_bv_select.load(in, &np_bv);
                        }

                        // create the new node-to-path indexes
                        index_node_to_path(temp_file::create());
                    }
                    else {
                        // we're in the more recent encoding, so we can load
                        // the node-to-path indexes directly

                        np_bv.load(in);
                        np_bv_select.load(in, &np_bv);
                        np_iv.load(in);
                        nr_iv.load(in);
                        nx_iv.load(in);
                    }

                    if (file_version >= 6 && file_version <= 10) {
                        // load and ignore the component path set indexes (which have
                        // now been exported)
                        {
                            sdsl::int_vector<> path_ranks_iv;
                            path_ranks_iv.load(in);
                        }
                        {
                            sdsl::bit_vector path_ranks_bv;
                            path_ranks_bv.load(in);
                        }
                    }

                    if (file_version <= 11) {
                        // load and ignore the gPBWT entity vectors
                        {
                            sdsl::vlc_vector<> h_civ;
                            h_civ.load(in);
                        }
                        // and the thread starts
                        {
                            sdsl::vlc_vector<> ts_civ;
                            ts_civ.load(in);
                        }
                        // and the B arrays
                        {
                            sdsl::wt_rlmn<sdsl::sd_vector<>> bs_single_array;
                            bs_single_array.load(in);
                        }
                    }
                }
                    break;
                default:
                    throw XGFormatError("Unimplemented XG format version: " + std::to_string(file_version));
            }
        } catch (const XGFormatError& e) {
            // Pass XGFormatErrors through
            throw e;
        } catch (const std::bad_alloc& e) {
            // We get std::bad_alloc generally if we try to read arbitrary data as an xg index.
            throw XGFormatError("XG input data not in XG version " + std::to_string(file_version) + " format (" + e.what() + ")");
        } catch (const std::exception& e) {
            // Other things will get re-thrown with a hint.
            std::cerr << "error [xg]: Unexpected error parsing XG data. Is it in version " << file_version << " XG format?" << std::endl;
            throw e;
        }
    }


    ////////////////////////////////////////////////////////////////////////////
    // Here is XPPath
    ////////////////////////////////////////////////////////////////////////////

    size_t XPPath::step_rank_at_position(size_t pos) const {
        return offsets_rank(pos+1)-1;
    }

    // TODO @ekg Do I need to write out all of these for our purpose?
    size_t XPPath::serialize(std::ostream& out,
                             sdsl::structure_tree_node* v,
                             std::string name) const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written = 0;
        written += sdsl::write_member(min_handle, out, child, "min_handle" + name);
        written += handles.serialize(out, child, "path_handles_" + name);
        written += offsets.serialize(out, child, "path_node_starts_" + name);
        written += offsets_rank.serialize(out, child, "path_node_starts_rank_" + name);
        written += offsets_select.serialize(out, child, "path_node_starts_select_" + name);
        written += sdsl::write_member(is_circular, out, child, "is_circular_" + name);
        sdsl::structure_tree::add_size(child, written);
        return written;
    }

    XPPath::XPPath(const std::string& path_name,
                   const std::vector<handle_t>& path,
                   bool is_circular,
                   XG& graph) {

#ifdef debug_path_index
        std::cerr << "Constructing xgpath for path with handles:" << std::endl;
    for (handle_t visiting : path) {
        std::cerr << "\t" << as_integer(visiting) << std::endl;
    }
#endif

        // The circularity flag is just a normal bool
        this->is_circular = is_circular;

        // handle integer values, the literal path
        sdsl::int_vector<> handles_iv;
        sdsl::util::assign(handles_iv, sdsl::int_vector<>(path.size()));
        // directions of traversal (typically forward, but we allow backwards)
        sdsl::bit_vector directions_bv;
        sdsl::util::assign(directions_bv, sdsl::bit_vector(path.size()));

        size_t path_off = 0;
        size_t members_off = 0;
        size_t positions_off = 0;
        size_t path_length = 0;
        // Start out with the max integer, as a handle, as our minimum-valued handle in the path.
        uint64_t min_handle_int = (path.size() ? as_integer(path[0]) : 0);

        // determine min handle value which occurs
        for (size_t i = 1; i < path.size(); ++i) {
            if (as_integer(path[i]) < min_handle_int) {
                min_handle_int = as_integer(path[i]);
            }
        }
        min_handle = as_handle(min_handle_int);

#ifdef debug_path_index
        std::cerr << "Basing on minimum handle value " << as_integer(min_handle) << " (aka " << min_handle_int << ")" << std::endl;
#endif

        // determine total length and record handles
        for (size_t i = 0; i < path.size(); ++i) {
            const handle_t& handle = path[i];
            path_length += graph.get_length(handle);
            handles_iv[i] = as_integer(local_handle(handle));

#ifdef debug_path_index
            std::cerr << "Recorded handle as " << handles_iv[i] << std::endl;
#endif

            // we will explode if the node isn't in the graph
        }
        sdsl::util::assign(handles, sdsl::enc_vector<>(handles_iv));

#ifdef debug_path_index
        for (size_t i = 0; i < path.size(); i++) {
        std::cerr << "Encoded handle as " << handles[i] << std::endl;
    }
#endif

        // make the bitvector for path offsets
        sdsl::bit_vector offsets_bv;
        sdsl::util::assign(offsets_bv, sdsl::bit_vector(path_length));

        //cerr << "path " << path_name << " has " << path.size() << endl;
        for (size_t i = 0; i < path.size(); ++i) {
            //cerr << i << endl;
            auto& handle = path[i];
            // record position of node
            offsets_bv[path_off] = 1;
            // and update the offset counter
            path_off += graph.get_length(handle);
        }
        // and path offsets
        sdsl::util::assign(offsets, sdsl::rrr_vector<>(offsets_bv));
        // and set up rank/select dictionary on them
        sdsl::util::assign(offsets_rank, sdsl::rrr_vector<>::rank_1_type(&offsets));
        sdsl::util::assign(offsets_select, sdsl::rrr_vector<>::select_1_type(&offsets));
    }

    void XPPath::load(std::istream& in) {
        sdsl::read_member(min_handle, in);
        handles.load(in);
        offsets.load(in);
        offsets_rank.load(in, &offsets);
        offsets_select.load(in, &offsets);
        sdsl::read_member(is_circular, in);
    }
}