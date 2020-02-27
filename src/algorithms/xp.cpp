#include "xp.hpp"

#include <arpa/inet.h>

namespace xp {

    using namespace handlegraph;

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

    void XP::from_enumerators(const std::function<void(const std::function<void(const std::string& seq, const nid_t& node_id)>&)>& for_each_sequence,
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

    }



    // FIXME This was directly copied from xg.cpp
    uint32_t XP::get_magic_number(void) const {
        return 4143290017ul;
    }

    void XP::serialize_members(ostream& out) const {
        serialize_and_measure(out);
    }

    // TODO Clean this up.
    /**
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
    **/

    void XP::deserialize_members(std::istream& in) {
        // simple alias to match an external interface
        load(in);
    }

    // FIXME This does not work as it is.
    /**
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
     **/


    ////////////////////////////////////////////////////////////////////////////
    // Here is XPPath
    ////////////////////////////////////////////////////////////////////////////

    size_t XPPath::step_rank_at_position(size_t pos) const {
        return offsets_rank(pos+1)-1;
    }

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

namespace temp_file {

// We use this to make the API thread-safe
    std::recursive_mutex monitor;

    std::string temp_dir;

/// Because the names are in a static object, we can delete them when
/// std::exit() is called.
    struct Handler {
        std::set<std::string> filenames;
        std::string parent_directory;
        ~Handler() {
            // No need to lock in static destructor
            for (auto& filename : filenames) {
                std::remove(filename.c_str());
            }
            if (!parent_directory.empty()) {
                // There may be extraneous files in the directory still (like .fai files)
                auto directory = opendir(parent_directory.c_str());

                dirent* dp;
                while ((dp = readdir(directory)) != nullptr) {
                    // For every item still in it, delete it.
                    // TODO: Maybe eventually recursively delete?
                    std::remove((parent_directory + "/" + dp->d_name).c_str());
                }
                closedir(directory);

                // Delete the directory itself
                std::remove(parent_directory.c_str());
            }
        }
    } handler;

    std::string create(const std::string& base) {
        std::lock_guard<recursive_mutex> lock(monitor);

        if (handler.parent_directory.empty()) {
            // Make a parent directory for our temp files
            string tmpdirname_cpp = get_dir() + "/xp-XXXXXX";
            char* tmpdirname = new char[tmpdirname_cpp.length() + 1];
            strcpy(tmpdirname, tmpdirname_cpp.c_str());
            auto got = mkdtemp(tmpdirname);
            if (got != nullptr) {
                // Save the directory we got
                handler.parent_directory = got;
            } else {
                cerr << "[xp]: couldn't create temp directory: " << tmpdirname << endl;
                exit(1);
            }
            delete [] tmpdirname;
        }

        std::string tmpname = handler.parent_directory + "/" + base + "XXXXXX";
        // hack to use mkstemp to get us a safe temporary file name
        int fd = mkstemp(&tmpname[0]);
        if(fd != -1) {
            // we don't leave it open; we are assumed to open it again externally
            close(fd);
        } else {
            cerr << "[xp]: couldn't create temp file on base "
                 << base << " : " << tmpname << endl;
            exit(1);
        }
        handler.filenames.insert(tmpname);
        return tmpname;
    }

    std::string create() {
        // No need to lock as we call this thing that locks
        return create("xp-");
    }

    void remove(const std::string& filename) {
        std::lock_guard<recursive_mutex> lock(monitor);

        std::remove(filename.c_str());
        handler.filenames.erase(filename);
    }

    void set_dir(const std::string& new_temp_dir) {
        std::lock_guard<recursive_mutex> lock(monitor);

        temp_dir = new_temp_dir;
    }

    std::string get_dir() {
        std::lock_guard<recursive_mutex> lock(monitor);

        // Get the default temp dir from environment variables.
        if (temp_dir.empty()) {
            const char* system_temp_dir = nullptr;
            for(const char* var_name : {"TMPDIR", "TMP", "TEMP", "TEMPDIR", "USERPROFILE"}) {
                if (system_temp_dir == nullptr) {
                    system_temp_dir = getenv(var_name);
                }
            }
            temp_dir = (system_temp_dir == nullptr ? "/tmp" : system_temp_dir);
        }

        return temp_dir;
    }

}