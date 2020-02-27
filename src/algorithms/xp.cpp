#include "xp.hpp"

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
    void XP::from_handle_graph(const PathHandleGraph& graph) {
        // the graph must be compacted for this to work
        std::vector<uint64_t> position_map(graph.get_node_count()+1);
        uint64_t len = 0;
        graph.for_each_handle([&](const handle_t& h) {
            position_map[number_bool_packing::unpack_number(h)] = len;
            uint64_t hl = graph.get_length(h);
            len += hl;
        });
        position_map[position_map.size()-1] = len;
        std::cout << "The current graph to index has nucleotide length: " << len << std::endl;

        graph.for_each_path_handle([&](const path_handle_t& path) {
            std::vector<handle_t> p;
            graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
                handle_t h = graph.get_handle_of_step(occ);
                p.push_back(h);
                uint64_t hl = graph.get_length(h);
            });
            std::string path_name = graph.get_path_name(path);
            std::cout << "Indexing: " << path_name << std::endl;
            XPPath* path_index = new XPPath(path_name, p, false, graph);
        });
    }

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
                   const handlegraph::PathHandleGraph& graph) {

#ifdef debug_path_index
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

    handle_t XPPath::local_handle(const handle_t& handle) const {
        if (as_integer(handle) < as_integer(min_handle)) {
            throw std::runtime_error("Handle with value " + std::to_string(as_integer(handle)) +
                                     " cannot be converted to local space based at min handle with value " + std::to_string(as_integer(min_handle)));
        } else {
            return as_handle(as_integer(handle)-as_integer(min_handle));
        }
    }
}