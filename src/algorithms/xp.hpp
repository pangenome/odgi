#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <dirent.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_arrays.hpp>

#include <handlegraph/util.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_position_handle_graph.hpp>

namespace xp {

    using nid_t = handlegraph::nid_t;

    class XPPath;

    /**
 * Thrown when attempting to interpret invalid data as an XP index.
 */
    class XPFormatError : public std::runtime_error {
        // Use the runtime_error constructor
        using std::runtime_error::runtime_error;
    };

    class XPQueryError : public std::runtime_error {
        // Use the runtime_error constructor
        using std::runtime_error::runtime_error;
    };

    /**
    * Provides succinct storage for the positional paths of a graph.
    */
    class XP {
    public:

        ////////////////////////////////////////////////////////////////////////////
        // Here are the ways we can construct an XP object from a graph
        ////////////////////////////////////////////////////////////////////////////

        XP(void) = default;

        ~XP(void);

        // We cannot move, assign, or copy until we add code to point sdsl supports
        // at the new addresses for their vectors.
        XP(const XP &other) = delete;

        XP(XP &&other) = delete;

        XP &operator=(const XP &other) = delete;

        XP &operator=(XP &&other) = delete;

        ////////////////////////////////////////////////////////////////////////////
        // Here is the handle graph API
        ////////////////////////////////////////////////////////////////////////////

        /// Build the path index from a simple graph.
        void from_handle_graph(const handlegraph::PathHandleGraph &graph);

        /// Load this XP index from a stream. Throw an XPFormatError if the stream
        /// does not produce a valid XP file.
        void load(std::istream &in);

        /// Alias for load() to match the SerializableHandleGraph interface.
        void deserialize_members(std::istream &in);

        /// Write this XP index to a stream.
        size_t
        serialize_and_measure(std::ostream &out, sdsl::structure_tree_node *s = nullptr, std::string name = "") const;

        /// Alias for serialize_and_measure().
        void serialize_members(std::ostream &out) const;

        /// Is this path in the index?
        bool has_path(const std::string& path_name) const;

        bool has_position(const std::string& path_name, size_t nuc_pos) const;

        /// Look up the path handle for the given path name
        handlegraph::path_handle_t get_path_handle(const std::string &path_name) const;

        /// Get a node handle (node ID and orientation) from a handle to a step on a path
        handlegraph::handle_t get_handle_of_step(const handlegraph::step_handle_t& step_handle) const;

        /// Returns the total length of sequence in the path
        size_t get_path_length(const handlegraph::path_handle_t& path_handle) const;

        /// Returns a handle to the path that an step is on
        handlegraph::path_handle_t get_path_handle_of_step(const handlegraph::step_handle_t& step_handle) const;

        /// Get the step at a given position
        handlegraph::step_handle_t
        get_step_at_position(const handlegraph::path_handle_t &path, const size_t &position) const;

        /// Look up the name of a path from a handle to it
        std::string get_path_name(const handlegraph::path_handle_t &path_handle) const;

        /// Look up the pangenome position by given path name and nucleotide position
        /// Returns 0 if the given path name is not in the index.
        /// Returns 0 if the given position is not in the given path.
        size_t get_pangenome_pos(const std::string &path_name, const size_t &nuc_pos) const;

        /// Get the path of the given path name
        const XPPath& get_path(const std::string& name) const;

        std::vector<XPPath *> get_paths() const;

        const sdsl::enc_vector<>& get_pos_map_iv() const;

        const sdsl::int_vector<>& get_pn_iv() const;

        size_t path_count = 0;

        char start_marker = '#';
        char end_marker = '$';

    private:
        ////////////////////////////////////////////////////////////////////////////
        // Here is path storage
        ////////////////////////////////////////////////////////////////////////////

        // path names
        sdsl::int_vector<> pn_iv; // path names
        sdsl::csa_wt<> pn_csa; // path name compressed suffix array
        sdsl::bit_vector pn_bv;  // path name starts in uncompressed version of csa
        sdsl::rank_support_v<1> pn_bv_rank;
        sdsl::bit_vector::select_1_type pn_bv_select;
        sdsl::int_vector<> pi_iv; // path ids by rank in the path names

        sdsl::enc_vector<> pos_map_iv; // store each offset of each node in the sequence vector

        std::vector<XPPath *> paths; // path structure
    };

    class XPPath {
    public:
        XPPath() = default;

        ~XPPath() = default;

        // Path name is required here only for complaining intelligently when
        // something goes wrong. We can also spit out the total unique members,
        // because in here is the most efficient place to count them.
        XPPath(const std::string &path_name,
               const std::vector<handlegraph::handle_t> &path,
               bool is_circular,
               const handlegraph::PathHandleGraph &graph);
        // Path names are stored in the XP object, in a compressed fashion, and are
        // not duplicated here.

        // These contain rank and select supports and so cannot move or be copied
        // without code to update them.
        XPPath(const XPPath &other) = delete;

        XPPath(XPPath &&other) = delete;

        XPPath &operator=(const XPPath &other) = delete;

        XPPath &operator=(XPPath &&other) = delete;

        handlegraph::handle_t min_handle;

        sdsl::enc_vector<> handles;
        sdsl::rrr_vector<> offsets;
        sdsl::rrr_vector<>::rank_1_type offsets_rank;
        sdsl::rrr_vector<>::select_1_type offsets_select;
        bool is_circular = false;

        void load(std::istream &in);

        size_t serialize(std::ostream &out,
                         sdsl::structure_tree_node *v = nullptr,
                         std::string name = "") const;

        size_t step_rank_at_position(size_t pos) const;

        handlegraph::handle_t local_handle(const handlegraph::handle_t &handle) const;
        handlegraph::handle_t external_handle(const handlegraph::handle_t& handle) const;
        handlegraph::handle_t handle(size_t offset) const;

    };

    /**
 * Temporary files. Create with create() and remove with remove(). All
 * temporary files will be deleted when the program exits normally or with
 * std::exit(). The files will be created in a directory determined from
 * environment variables, though this can be overridden with set_dir().
 * The interface is thread-safe.
 */
    namespace temp_file {

/// Create a temporary file starting with the given base name
        std::string create(const std::string &base);

/// Create a temporary file
        std::string create();

/// Remove a temporary file
        void remove(const std::string &filename);

/// Set a temp dir, overriding system defaults and environment variables.
        void set_dir(const std::string &new_temp_dir);

/// Get the current temp dir
        std::string get_dir();

    } // namespace temp_file

}
