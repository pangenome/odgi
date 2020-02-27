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
 * Thrown when attempting to interpret invalid data as an XP index.
 */
    class XPFormatError : public std::runtime_error {
        // Use the runtime_error constructor
        using std::runtime_error::runtime_error;
    };
 /**
 * Provides succinct storage for the positional paths of a graph.
 */
    // class XP : public handlegraph::PathPositionHandleGraph, public handlegraph::SerializableHandleGraph, public handlegraph::VectorizableHandleGraph {
    // TODO @ekg Do I need any of these abstract classes? I would have to implement a huge amount of functions. See path_handle_graph.hpp for example.
    class XP {
    public:

        ////////////////////////////////////////////////////////////////////////////
        // Here are the ways we can construct an XP object from a graph
        ////////////////////////////////////////////////////////////////////////////

        XP() = default;

        ~XP();

        // We cannot move, assign, or copy until we add code to point sdsl supports
        // at the new addresses for their vectors.
        XP(const XP &other) = delete;

        XP(XP&& other) = delete;

        XP &operator=(const XP &other) = delete;

        XP &operator=(XP &&other) = delete;

        // General public statistics
        size_t seq_length = 0;
        size_t node_count = 0;
        size_t edge_count = 0;
        size_t path_count = 0;

        ////////////////////////////////////////////////////////////////////////////
        // Here is the handle graph API
        ////////////////////////////////////////////////////////////////////////////

        /// Look up the handle for the node with the given ID in the given orientation
        virtual handlegraph::handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const;

        size_t id_to_rank(const nid_t& id) const;

        // Build the path index from a simple graph.
        void from_handle_graph(const handlegraph::HandleGraph &graph);

        /// Use external enumerators to drive graph construction.
        /// The order in which nodes are enumerated becomes the XP's node order.
        /// Note that we will get the best efficiency if the graph is enumerated in topological order.
        void from_enumerators(const std::function<void(const std::function<void(const std::string& seq, const nid_t& node_id)>&)>& for_each_sequence,
                              const std::function<void(const std::function<void(const nid_t& from, const bool& from_rev,
                                                                                const nid_t& to, const bool& to_rev)>&)>& for_each_edge,
                              const std::function<void(const std::function<void(const std::string& path_name,
                                                                                const nid_t& node_id, const bool& is_rev,
                                                                                const std::string& cigar, const bool& is_empty,
                                                                                const bool& is_circular)>&)>& for_each_path_element,
                              bool validate = false, std::string basename = "");

        // Get our magic number
        uint32_t get_magic_number() const;

        // Load this XP index from a stream. Throw an XPFormatError if the stream
        // does not produce a valid XP file.
        void load(std::istream& in);

        // Alias for load() to match the SerializableHandleGraph interface.
        void deserialize_members(std::istream& in);

        // Write this XP index to a stream.
        size_t serialize_and_measure(std::ostream& out, sdsl::structure_tree_node* s = nullptr, std::string name = "") const;
        // Alias for serialize_and_measure().
        void serialize_members(std::ostream& out) const;

        /// Look up the path handle for the given path name
        handlegraph::path_handle_t get_path_handle(const std::string& path_name) const;

        /// Get the step at a given position
        handlegraph::step_handle_t get_step_at_position(const handlegraph::path_handle_t& path, const size_t& position) const;

        /// Get the length of a node
        virtual size_t get_length(const handlegraph::handle_t& handle) const;

        char start_marker = '#';
        char end_marker = '$';

    private:
        ////////////////////////////////////////////////////////////////////////////
        // Here is the New Way (locally traversable graph storage)
        // Everything should be rewritten in terms of these members
        ////////////////////////////////////////////////////////////////////////////

        /// locally traversable graph storage
        ///
        /// Encoding designed for efficient compression, cache locality, and relativistic traversal of the graph.
        ///
        /// node := { header, edges_to, edges_from }
        /// header := { node_id, node_start, node_length, edges_to_count, edges_from_count }
        /// node_id := integer
        /// node_start := integer (offset in s_iv)
        /// node_length := integer
        /// edges_to_count := integer
        /// edges_from_count := integer
        /// edges_to := { edge_to, ... }
        /// edges_from := { edge_from, ... }
        /// edge_to := { offset_to_previous_node, edge_type }
        /// edge_to := { offset_to_next_node, edge_type }
        sdsl::int_vector<> g_iv;
        /// delimit node records to allow lookup of nodes in g_civ by rank
        sdsl::bit_vector g_bv;
        sdsl::rank_support_v<1> g_bv_rank;
        sdsl::bit_vector::select_1_type g_bv_select;

        // Let's define some offset ints
        const static int G_NODE_ID_OFFSET = 0;
        const static int G_NODE_SEQ_START_OFFSET = 1;
        const static int G_NODE_LENGTH_OFFSET = 2;
        const static int G_NODE_TO_COUNT_OFFSET = 3;
        const static int G_NODE_FROM_COUNT_OFFSET = 4;
        const static int G_NODE_HEADER_LENGTH = 5;

        const static int G_EDGE_OFFSET_OFFSET = 0;
        const static int G_EDGE_TYPE_OFFSET = 1;
        const static int G_EDGE_LENGTH = 2;

        // And the edge types (so we don't confuse our magic numbers)
        const static int EDGE_TYPE_MIN = 1;
        const static int EDGE_TYPE_END_START = 1;
        const static int EDGE_TYPE_END_END = 2;
        const static int EDGE_TYPE_START_START = 3;
        const static int EDGE_TYPE_START_END = 4;
        const static int EDGE_TYPE_MAX = 4;

        ////////////////////////////////////////////////////////////////////////////
        // And here are the bits for tracking actual node IDs
        ////////////////////////////////////////////////////////////////////////////

        // maintain old ids from input, ranked as in s_iv and s_bv
        int64_t min_id = 0; // id ranges don't have to start at 0
        int64_t max_id = 0;
        sdsl::int_vector<> r_iv; // ids-id_min is the rank

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

        // TODO @ekg I won't need any of this, right?
        /**
        // node->path membership
        sdsl::int_vector<> np_iv;
        // node->path rank
        sdsl::int_vector<> nr_iv;
        // node->path position/orientation
        sdsl::int_vector<> nx_iv;
        sdsl::bit_vector np_bv; // entity delimiters in both vectors
        //sdsl::rank_support_v<1> np_bv_rank;
        sdsl::bit_vector::select_1_type np_bv_select;
        **/
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
               XP& graph);
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
        //sdsl::rrr_vector directions; // forward or backward through nodes
        sdsl::rrr_vector<> offsets;
        sdsl::rrr_vector<>::rank_1_type offsets_rank;
        sdsl::rrr_vector<>::select_1_type offsets_select;
        bool is_circular = false;

        void load(std::istream& in);

        size_t serialize(std::ostream& out,
                         sdsl::structure_tree_node* v = NULL,
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
