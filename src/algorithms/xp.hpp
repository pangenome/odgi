#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <queue>
#include <omp.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <functional>
#include <utility>
#include <tuple>
#include <sys/types.h>
#include <dirent.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/dac_vector.hpp>
#include <sdsl/vlc_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_arrays.hpp>

#include <handlegraph/types.hpp>
#include <handlegraph/iteratee.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_position_handle_graph.hpp>
#include <handlegraph/serializable_handle_graph.hpp>

namespace xp {

    using nid_t = handlegraph::nid_t;

    class XPPath;
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

        XP(XP&& other) = delete;

        XP &operator=(const XP &other) = delete;

        XP &operator=(XP &&other) = delete;

        ////////////////////////////////////////////////////////////////////////////
        // Here is the handle graph API
        ////////////////////////////////////////////////////////////////////////////

        /// Build the path index from a simple graph.
        void from_handle_graph(const handlegraph::PathHandleGraph& graph);

        /// Look up the handle for the node with the given ID in the given orientation
        // TODO If not implemented, the linker crashes because of virtual declaration.
        // virtual handlegraph::handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const;

        size_t id_to_rank(const nid_t& id) const;

        /// Load this XP index from a stream. Throw an XPFormatError if the stream
        /// does not produce a valid XP file.
        void load(std::istream& in);

        /// Alias for load() to match the SerializableHandleGraph interface.
        void deserialize_members(std::istream& in);

        /// Write this XP index to a stream.
        size_t serialize_and_measure(std::ostream& out, sdsl::structure_tree_node* s = nullptr, std::string name = "") const;
        /// Alias for serialize_and_measure().
        void serialize_members(std::ostream& out) const;

        /// Look up the path handle for the given path name
        handlegraph::path_handle_t get_path_handle(const std::string& path_name) const;

        /// Get the step at a given position
        handlegraph::step_handle_t get_step_at_position(const handlegraph::path_handle_t& path, const size_t& position) const;

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

        sdsl::int_vector<> position_map; // store each offset of each node in the sequence vector

        std::vector<XPPath*> paths; // path structure
    };

    class XPPath {
    public:
        XPPath() = default;
        ~XPPath() = default;
        // Path name is required here only for complaining intelligently when
        // something goes wrong. We can also spit out the total unique members,
        // because in here is the most efficient place to count them.
        XPPath(const std::string& path_name,
               const std::vector<handlegraph::handle_t>& path,
               bool is_circular,
               const handlegraph::PathHandleGraph& graph);
        // Path names are stored in the XP object, in a compressed fashion, and are
        // not duplicated here.

        // These contain rank and select supports and so cannot move or be copied
        // without code to update them.
        XPPath(const XPPath& other) = delete;
        XPPath(XPPath&& other) = delete;
        XPPath& operator=(const XPPath& other) = delete;
        XPPath& operator=(XPPath&& other) = delete;
        handlegraph::handle_t min_handle;

        sdsl::enc_vector<> handles;
        sdsl::rrr_vector<> offsets;
        sdsl::rrr_vector<>::rank_1_type offsets_rank;
        sdsl::rrr_vector<>::select_1_type offsets_select;
        bool is_circular = false;

        void load(std::istream& in);

        size_t serialize(std::ostream& out,
                         sdsl::structure_tree_node* v = nullptr,
                         std::string name = "") const;

        size_t step_rank_at_position(size_t pos) const;

        handlegraph::handle_t local_handle(const handlegraph::handle_t& handle) const;
    };

}

/**
 * Temporary files. Create with create() and remove with remove(). All
 * temporary files will be deleted when the program exits normally or with
 * std::exit(). The files will be created in a directory determined from
 * environment variables, though this can be overridden with set_dir().
 * The interface is thread-safe.
 */
namespace temp_file {

/// Create a temporary file starting with the given base name
    std::string create(const std::string& base);

/// Create a temporary file
    std::string create();

/// Remove a temporary file
    void remove(const std::string& filename);

/// Set a temp dir, overriding system defaults and environment variables.
    void set_dir(const std::string& new_temp_dir);

/// Get the current temp dir
    std::string get_dir();

} // namespace temp_file
