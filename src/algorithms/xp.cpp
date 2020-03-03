#include "xp.hpp"

#include <arpa/inet.h>
#include <mutex>

namespace xp {

    using namespace handlegraph;

    ////////////////////////////////////////////////////////////////////////////
    // Here is XP
    ////////////////////////////////////////////////////////////////////////////

    XP::~XP(void) {
        // Clean up any created XPPaths
        while (!paths.empty()) {
            delete paths.back();
            paths.pop_back();
        }
    }

    /// build the graph from a graph handle
    void XP::from_handle_graph(const PathHandleGraph &graph) {
        // create temporary file for path names
        std::string basename;
        if (basename.empty()) {
            basename = temp_file::create();
        }
        std::string path_names;
        // the graph must be compacted for this to work
        sdsl::int_vector<> position_map;
        sdsl::util::assign(position_map, sdsl::int_vector<>(graph.get_node_count() + 1));
        uint64_t len = 0;
        graph.for_each_handle([&](const handle_t &h) {
            position_map[number_bool_packing::unpack_number(h)] = len;
            uint64_t hl = graph.get_length(h);
            len += hl;
        });
        position_map[position_map.size() - 1] = len;
#ifdef debug_from_handle_graph
        std::err << "[XP CONSTRUCTION]: The current graph to index has nucleotide length: " << len << std::endl;
        std::err << "[XP CONSTRUCTION]: position_map: ";
        for (size_t i = 0; i < position_map.size() - 1; i++) {
            std::cerr << position_map[i] << ",";
        }
        std::err << position_map[position_map.size() - 1] << std::endl;
#endif

        graph.for_each_path_handle([&](const path_handle_t &path) {
            std::vector<handle_t> p;
            graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                handle_t h = graph.get_handle_of_step(occ);
                p.push_back(h);
                uint64_t hl = graph.get_length(h);
            });
            std::string path_name = graph.get_path_name(path);
            std::cout << "[XP CONSTRUCTION]: Indexing path: " << path_name << std::endl;
            XPPath *path_index = new XPPath(path_name, p, false, graph);
            // Add path_index to paths.
            paths.push_back(path_index);
            path_names += start_marker + path_name + end_marker;
        });
        // assign the position map iv
        sdsl::util::assign(pos_map_iv, sdsl::enc_vector<>(position_map));
        // set the path counts
        path_count = paths.size();
#ifdef from_handle_graph
        std::cout << "[XP CONSTRUCTION]: path_count: " << path_count << std::endl;
#endif
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

        // write path names to temp file
        std::string path_name_file = basename + ".pathnames.iv";
        sdsl::store_to_file((const char *) path_names.c_str(), path_name_file);
        // read file and construct compressed suffix array
        sdsl::construct(pn_csa, path_name_file, 1);
    }

    std::string XP::get_path_name(const handlegraph::path_handle_t &path_handle) const {
        uint64_t rank = as_integer(path_handle);
        size_t start = pn_bv_select(rank) + 1; // step past '#'
        size_t end = rank == path_count ? pn_iv.size() : pn_bv_select(rank + 1);
        end -= 1;  // step before '$'
        std::string name;
        name.resize(end - start);
        for (size_t i = start; i < end; ++i) {
            name[i - start] = pn_iv[i];
        }
        return name;
    }

    void XP::serialize_members(std::ostream &out) const {
        serialize_and_measure(out);
    }

    size_t XP::serialize_and_measure(std::ostream &out, sdsl::structure_tree_node *s, std::string name) const {

        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
        size_t written = 0;

        // Do the magic number
        out << "XP";
        written += 2;

        // POSITION MAP STUFF
        written += pos_map_iv.serialize(out, child, "position_map");
        // PATH STUFF
        written += sdsl::write_member(path_count, out, child, "path_count");

        // Treat the paths as their own node
        size_t paths_written = 0;
        auto paths_child = sdsl::structure_tree::add_child(child, "paths", sdsl::util::class_name(*this));

        paths_written += pn_iv.serialize(out, paths_child, "path_names");
        paths_written += pn_csa.serialize(out, paths_child, "path_names_csa");
        paths_written += pn_bv.serialize(out, paths_child, "path_names_starts");
        paths_written += pn_bv_rank.serialize(out, paths_child, "path_names_starts_rank");
        paths_written += pn_bv_select.serialize(out, paths_child, "path_names_starts_select");
        paths_written += pi_iv.serialize(out, paths_child, "path_ids");

        for (size_t i = 0; i < paths.size(); i++) {
            XPPath *path = paths[i];
            paths_written += path->serialize(out, paths_child,
                                             "path:" + XP::get_path_name(handlegraph::as_path_handle(i + 1)));
        }

        sdsl::structure_tree::add_size(paths_child, paths_written);
        written += paths_written;

        sdsl::structure_tree::add_size(child, written);
        return written;
    }

    void XP::deserialize_members(std::istream &in) {
        // simple alias to match an external interface
        load(in);
    }

    void XP::load(std::istream &in) {

        if (!in.good()) {
            throw XPFormatError("Index file does not exist or index stream cannot be read");
        }

        // We need to look for the magic value
        char buffer;
        in.get(buffer);
        if (buffer == 'X') {
            in.get(buffer);
            if (buffer == 'P') {
                // We found the magic value!

            } else {
                // Put back both characters
                in.unget();
                in.unget();
            }
        } else {
            // Put back the one character
            in.unget();
        }

        try {
            pos_map_iv.load(in);
            sdsl::read_member(path_count, in);
            pn_iv.load(in);
            pn_csa.load(in);
            pn_bv.load(in);
            pn_bv_rank.load(in, &pn_bv);
            pn_bv_select.load(in, &pn_bv);
            pi_iv.load(in);

            for (size_t i = 0; i < path_count; ++i) {
                auto path = new XPPath;
                // Load the path, giving it the file version and a
                // rank-to-ID comversion function for format upgrade
                // purposes.
                path->load(in);
                paths.push_back(path);
            }
        } catch (const XPFormatError &e) {
            // Pass XGFormatErrors through
            throw e;
        } catch (const std::bad_alloc &e) {
            // We get std::bad_alloc generally if we try to read arbitrary data as an xg index.
            throw XPFormatError("XP input data not in XP version " + std::to_string(42) + " format (" + e.what() + ")");
        } catch (const std::exception &e) {
            // Other things will get re-thrown with a hint.
            std::cerr << "error [xp]: Unexpected error parsing XP data. Is it in version " << "EASTER EGG"
                      << " XP format?" << std::endl;
            throw e;
        }
    }

    path_handle_t XP::get_path_handle(const std::string& path_name) const {
        // find the name in the csa
        std::string query = start_marker + path_name + end_marker;
        auto occs = locate(pn_csa, query);
        if (occs.size() > 1) {
            std::cerr << "error [xp]: multiple hits for " << query << std::endl;
            exit(1);
        }
        if(occs.size() == 0) {
            // This path does not exist. Give back 0, which can never be a real path
            // rank.
            return as_path_handle(0);
        }
        return as_path_handle(pn_bv_rank(occs[0])+1); // step past '#'
    }

    size_t XP::get_path_length(const path_handle_t& path_handle) const {
        return paths[as_integer(path_handle) - 1]->offsets.size();
    }

    /// Get the step at a given position
    step_handle_t XP::get_step_at_position(const path_handle_t& path, const size_t& position) const {
        if (position >= get_path_length(path)) {
            throw std::runtime_error("Cannot get position " + std::to_string(position) + " along path " +
                                get_path_name(path) + " of length " + std::to_string(get_path_length(path)));
        }

        const auto& xppath = *paths[as_integer(path) - 1];
        step_handle_t step;
        as_integers(step)[0] = as_integer(path);
        as_integers(step)[1] = xppath.step_rank_at_position(position);
        return step;
    }

    path_handle_t XP::get_path_handle_of_step(const step_handle_t& step_handle) const {
        return as_path_handle(as_integers(step_handle)[0]);
    }

    handle_t XP::get_handle_of_step(const step_handle_t& step_handle) const {
        const auto& xppath = *paths[as_integer(get_path_handle_of_step(step_handle)) - 1];
        return xppath.handle(as_integers(step_handle)[1]);
    }

    const XPPath& XP::get_path(const std::string &name) const {
        handlegraph::path_handle_t p_h = get_path_handle(name);
        return *paths[as_integer(p_h) - 1];
    }

    size_t XP::get_pangenome_pos(const std::string &path_name, const size_t &nuc_pos) const {
#ifdef debug_get_pangenome_pos
        std::cerr << "[GET_PANGENOME_POS]: path_name: " << path_name << std::endl;
#endif
        const XPPath& xppath = get_path(path_name);
        // TODO Is the nucleotide position there?!
        if (xppath.offsets.size() < nuc_pos) {
            return 0;
        }
        size_t step_rank = xppath.step_rank_at_position(nuc_pos - 1);
#ifdef debug_get_pangenome_pos
        std::cerr << "[GET_PANGENOME_POS]: step_rank: " << step_rank << std::endl;
#endif
        // get the path handle of the given path name
        handlegraph::path_handle_t p_h = get_path_handle(path_name);
#ifdef debug_get_pangenome_pos
        std::cerr << "[GET_PANGENOME_POS]: p_h: " << as_integer(p_h) << std::endl;
#endif
        // TODO Is the given path name even in the index?! Move this more to the top so it breaks before checking the nuc_pos.
        if (p_h == as_path_handle(0)) {
            std::cerr << "The given path name " << path_name << " is not in the index." << std::endl;
            exit(1);
        }
        step_handle_t step_pos = get_step_at_position(p_h, nuc_pos - 1);
#ifdef debug_get_pangenome_pos
        std::cerr << "[GET_PANGENOME_POS]: step_pos: " << as_integers(step_pos)[0] << as_integers(step_pos)[1] << std::endl;
#endif
        handle_t p = get_handle_of_step(step_pos);
#ifdef debug_get_pangenome_pos
        std::cerr << "[GET_PANGENOME_POS]: handle_t of step_pos: " << as_integer(p) << std::endl;
#endif
        // p = as_handle(as_integer(p) + 1);

        // handle position
        uint64_t handle_pos = pos_map_iv[number_bool_packing::unpack_number(p)];
#ifdef debug_get_pangenome_pos
        std::cerr << "[GET_PANGENOME_POS]: handle_pos: " << handle_pos << std::endl;
#endif
        // length of the handle
        uint64_t next_handle_pos = pos_map_iv[number_bool_packing::unpack_number(p) + 1];
#ifdef debug_get_pangenome_pos
        std::cerr << "[GET_PANGENOME_POS]: next_handle_pos: " << next_handle_pos << std::endl;
#endif
        uint64_t node_length = next_handle_pos - handle_pos;
#ifdef debug_get_pangenome_pos
        std::cerr << "[GET_PANGENOME_POS]: node_length: " << node_length << std::endl;
#endif
        bool is_rev = number_bool_packing::unpack_bit(p);
        // Adjust this for both strands!!!
        uint64_t offset_in_handle = nuc_pos - as_integer(p);
#ifdef debug_get_pangenome_pos
        std::cerr << "[GET_PANGENOME_POS]: offset_in_handle: " << offset_in_handle << std::endl;
#endif
        if (is_rev) {
            offset_in_handle = node_length - offset_in_handle - 1;
        }
        size_t pos_in_pangenome = handle_pos + offset_in_handle;

        return pos_in_pangenome;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Here is XPPath
    ////////////////////////////////////////////////////////////////////////////

    size_t XPPath::step_rank_at_position(size_t pos) const {
        return offsets_rank(pos + 1) - 1;
    }

    size_t XPPath::serialize(std::ostream &out,
                             sdsl::structure_tree_node *v,
                             std::string name) const {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
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

    XPPath::XPPath(const std::string &path_name,
                   const std::vector<handle_t> &path,
                   bool is_circular,
                   const handlegraph::PathHandleGraph &graph) {

#ifdef debug_xppath
        std::cerr << "Constructing xppath for path with handles:" << std::endl;
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

#ifdef debug_xppath
        std::cerr << "Basing on minimum handle value " << as_integer(min_handle) << " (aka " << min_handle_int << ")" << std::endl;
#endif
        // determine total length and record handles
        for (size_t i = 0; i < path.size(); ++i) {
            const handle_t &handle = path[i];
            path_length += graph.get_length(handle);
            handles_iv[i] = as_integer(local_handle(handle));

#ifdef debug_xppath
            std::cerr << "Recorded handle as " << handles_iv[i] << std::endl;
#endif
            // we will explode if the node isn't in the graph
        }
        sdsl::util::assign(handles, sdsl::enc_vector<>(handles_iv));

#ifdef debug_xppath
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
            auto &handle = path[i];
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

#ifdef debug_xppath
        std::cerr << "[XPPATH CONSTRUCTION]: offsets_bv: ";
        for (size_t i = 0; i < offsets_bv.size(); i++) {
            std::cerr << offsets_bv[i];
        }
        std::cerr << std::endl;

        std::cerr << "[XPPATH CONSTRUCTION]: offsets: ";
        for (size_t i = 0; i < offsets.size(); i++) {
            std::cerr << offsets[i];
        }
        std::cerr << std::endl;

        std::cerr << "[XPPATH CONSTRUCTION]: offsets_rank: ";
        for (size_t i = 0; i < offsets_rank.size(); i++) {
            std::cerr << offsets_rank(i+1) - 1;
        }
        std::cerr << std::endl;

        std::cerr << "[XPPATH CONSTRUCTION]: offsets_select: ";
        for (size_t i = 1; i <= path.size(); ++i) {
            std::cerr << offsets_select(i) << ",";
        }
        std::cerr << std::endl;
#endif
    }

    void XPPath::load(std::istream &in) {
        sdsl::read_member(min_handle, in);
        handles.load(in);
        offsets.load(in);
        offsets_rank.load(in, &offsets);
        offsets_select.load(in, &offsets);
        sdsl::read_member(is_circular, in);

#ifdef debug_load_xppath
        std::cout << "[XPPATH LOAD]: offsets: ";
        for (size_t i = 0; i < offsets.size(); i++) {
            std::cout << offsets[i];
        }
        std::cout << std::endl;

        std::cout << "[XPPATH LOAD]: offsets_rank: ";
        for (size_t i = 0; i < offsets_rank.size(); i++) {
            std::cout << offsets_rank(i+1) - 1;
        }
        std::cout << std::endl;

        std::cout << "[XPPATH LOAD]: offsets_select: ";
        for (size_t i = 1; i <= offsets_select.size(); ++i) {
            std::cout << offsets_select(i) << ",";
        }
        std::cout << std::endl;
#endif
    }

    handle_t XPPath::local_handle(const handle_t &handle) const {
        if (as_integer(handle) < as_integer(min_handle)) {
            throw std::runtime_error("Handle with value " + std::to_string(as_integer(handle)) +
                                     " cannot be converted to local space based at min handle with value " +
                                     std::to_string(as_integer(min_handle)));
        } else {
            return as_handle(as_integer(handle) - as_integer(min_handle));
        }
    }

    handle_t XPPath::handle(size_t offset) const {
        return external_handle(as_handle(handles[offset]));
    }

    handle_t XPPath::external_handle(const handle_t& handle) const {
        return as_handle(as_integer(handle)+as_integer(min_handle));
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
                for (auto &filename : filenames) {
                    std::remove(filename.c_str());
                }
                if (!parent_directory.empty()) {
                    // There may be extraneous files in the directory still (like .fai files)
                    auto directory = opendir(parent_directory.c_str());

                    dirent *dp;
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

        std::string create(const std::string &base) {
            std::lock_guard<std::recursive_mutex> lock(monitor);

            if (handler.parent_directory.empty()) {
                // Make a parent directory for our temp files
                std::string tmpdirname_cpp = get_dir() + "/xp-XXXXXX";
                char *tmpdirname = new char[tmpdirname_cpp.length() + 1];
                std::strcpy(tmpdirname, tmpdirname_cpp.c_str());
                auto got = mkdtemp(tmpdirname);
                if (got != nullptr) {
                    // Save the directory we got
                    handler.parent_directory = got;
                } else {
                    std::cerr << "[xp]: couldn't create temp directory: " << tmpdirname << std::endl;
                    exit(1);
                }
                delete[] tmpdirname;
            }

            std::string tmpname = handler.parent_directory + "/" + base + "XXXXXX";
            // hack to use mkstemp to get us a safe temporary file name
            int fd = mkstemp(&tmpname[0]);
            if (fd != -1) {
                // we don't leave it open; we are assumed to open it again externally
                close(fd);
            } else {
                std::cerr << "[xp]: couldn't create temp file on base "
                          << base << " : " << tmpname << std::endl;
                exit(1);
            }
            handler.filenames.insert(tmpname);
            return tmpname;
        }

        std::string create() {
            // No need to lock as we call this thing that locks
            return create("xp-");
        }

        void remove(const std::string &filename) {
            std::lock_guard<std::recursive_mutex> lock(monitor);

            std::remove(filename.c_str());
            handler.filenames.erase(filename);
        }

        void set_dir(const std::string &new_temp_dir) {
            std::lock_guard<std::recursive_mutex> lock(monitor);

            temp_dir = new_temp_dir;
        }

        std::string get_dir() {
            std::lock_guard<std::recursive_mutex> lock(monitor);

            // Get the default temp dir from environment variables.
            if (temp_dir.empty()) {
                const char *system_temp_dir = nullptr;
                for (const char *var_name : {"TMPDIR", "TMP", "TEMP", "TEMPDIR", "USERPROFILE"}) {
                    if (system_temp_dir == nullptr) {
                        system_temp_dir = getenv(var_name);
                    }
                }
                temp_dir = (system_temp_dir == nullptr ? "/tmp" : system_temp_dir);
            }

            return temp_dir;
        }
    }
}
