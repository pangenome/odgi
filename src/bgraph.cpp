#include "bgraph.hpp"

namespace betagraph{
    BGraph::BGraph(){

    };
    BGraph::~BGraph(){

    };
    handle_t BGraph::get_handle(const handlegraph::nid_t& node_id, bool is_reverse) const{
        return hhelper.pack(node_id, is_reverse);
    }
    

    handlegraph::nid_t BGraph::get_id(const handle_t& handle) const{
        return hhelper.unpack_number(handle);
    }
    

    bool BGraph::get_is_reverse(const handle_t& handle) const{
        return hhelper.unpack_bit(handle);
    }
    

    handle_t BGraph::flip(const handle_t& handle) const{
        return hhelper.toggle_bit(handle);
    }
    

    size_t BGraph::get_length(const handle_t& handle) const{
        return graph.backer.at(get_id(handle)).length;
    }

    std::string BGraph::get_sequence(const handle_t& handle) const{
        return graph.backer.at(get_id(handle)).sequence;
    }
    
    // bool follow_edges(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const{

    // }
    
    void BGraph::for_each_handle(const std::function<bool(const handle_t&)>& iteratee, bool parallel) const{
        std::atomic<bool> keep_going(true);

        while(keep_going){

            //if (parallel){

            //}
            //else{
                for (auto g : graph.backer){
                    if(!iteratee(get_handle(g.first, g.second.orientation))){
                        keep_going = false;
                    }

                    if (!keep_going){
                        break;
                    }
                }
            //}
            
        }
    }
    
    size_t BGraph::node_size() const{
        return graph.backer.size();
    }
    
    handlegraph::nid_t BGraph::min_node_id() const{
        return graph.min_node_id;
    }

    handlegraph::nid_t BGraph::max_node_id() const{
        return graph.max_node_id;
    }
  
//     /// Get the number of edges on the right (go_left = false) or left (go_left
//     /// = true) side of the given handle. The default implementation is O(n) in
//     /// the number of edges returned, but graph implementations that track this
//     /// information more efficiently can override this method.
     size_t BGraph::get_degree(const handle_t& handle, bool go_left) const{
         if (go_left){
             return graph.backer.at(as_integer(handle)).in_edges.size();
         }
         else{
             return graph.backer.at(as_integer(handle)).out_edges.size();
         }
     }
  
    /// Get the locally forward version of a handle
    handle_t BGraph::forward(const handle_t& handle) const{
        handle_t h(handle);
        if (hhelper.unpack_bit(h)){
            hhelper.toggle_bit(h);
        }
        return h;
    }
    
    /// A pair of handles can be used as an edge. When so used, the handles have a
    /// canonical order and orientation.
    // edge_t BGraph::edge_handle(const handle_t& left, const handle_t& right) const{

    // }
    
//     /// Such a pair can be viewed from either inward end handle and produce the
//     /// outward handle you would arrive at.
//     handle_t traverse_edge_handle(const edge_t& edge, const handle_t& left) const;
    
/**
 * This is the interface for a handle graph that stores embedded paths.
 */
    
    ////////////////////////////////////////////////////////////////////////////
    // Path handle interface that needs to be implemented
    ////////////////////////////////////////////////////////////////////////////
    
    /// Determine if a path name exists and is legal to get a path handle for.
     bool BGraph::has_path(const std::string& path_name) const{
         return paths.has_path(path_name);
     }
    
    /// Look up the path handle for the given path name.
    /// The path with that name must exist.
    path_handle_t BGraph::get_path_handle(const std::string& path_name) const{
        return as_path_handle(paths.get_id(path_name));
    }
    
    /// Look up the name of a path from a handle to it
    std::string BGraph::get_path_name(const path_handle_t& path_handle) const{
        return paths.get_name(path_handle);
    }
    
    /// Returns the number of node steps in the path
    size_t BGraph::get_step_count(const path_handle_t& path_handle) const{
        return paths.get_path_step_count(path_handle);
    }

    /// Returns the number of paths stored in the graph
    size_t BGraph::get_path_count() const{
        return paths.paths.size();
    }
    
    /// Execute a function on each path in the graph
    // TODO: allow stopping early?
    void BGraph::for_each_path_handle(const std::function<void(const path_handle_t&)>& iteratee) const{

        std::cerr << "Not implemented" << std::endl;
        exit(1);
    }
    
    /// Get a node handle (node ID and orientation) from a handle to an step on a path
    // handle_t BGraph::get_step(const step_handle_t& step_handle) const{

    // }
    
//     /// Get a handle to the first step in a path.
//     /// The path MUST be nonempty.
//     step_handle_t get_first_step(const path_handle_t& path_handle) const;
    
//     /// Get a handle to the last step in a path
//     /// The path MUST be nonempty.
//     step_handle_t get_last_step(const path_handle_t& path_handle) const;
    
//     /// Returns true if the step is not the last occurence on the path, else false
//     bool has_next_step(const step_handle_t& step_handle) const;
    
//     /// Returns true if the step is not the first occurence on the path, else false
//     bool has_previous_step(const step_handle_t& step_handle) const;
    
//     /// Returns a handle to the next step on the path
//     step_handle_t get_next_step(const step_handle_t& step_handle) const;
    
//     /// Returns a handle to the previous step on the path
//     step_handle_t get_previous_step(const step_handle_t& step_handle) const;
    
//     /// Returns a handle to the path that an step is on
//     path_handle_t get_path_handle_of_step(const step_handle_t& step_handle) const;
    
//     /// Returns the 0-based ordinal rank of a step on a path
//     size_t get_ordinal_rank_of_step(const step_handle_t& step_handle) const;

//     ////////////////////////////////////////////////////////////////////////////
//     // Additional optional interface with a default implementation
//     ////////////////////////////////////////////////////////////////////////////

//     /// Returns true if the given path is empty, and false otherwise
//     bool is_empty(const path_handle_t& path_handle) const;

//     ////////////////////////////////////////////////////////////////////////////
//     // Concrete utility methods
//     ////////////////////////////////////////////////////////////////////////////

//     /// Loop over all the steps along a path, from first through last
//     void for_each_step_in_path(const path_handle_t& path, const std::function<void(const step_handle_t&)>& iteratee) const;

// /**
//  * This is the interface for a handle graph that supports modification.
//  */
//     /*
//      * Note: All operations may invalidate path handles and step handles.
//      */
    
//     /// Create a new node with the given sequence and return the handle.
//     handle_t create_handle(const std::string& sequence);

//     /// Create a new node with the given id and sequence, then return the handle.
//     handle_t create_handle(const std::string& sequence, const handlegraph::nid_t& id);
    
//     /// Remove the node belonging to the given handle and all of its edges.
//     /// Does not update any stored paths.
//     /// Invalidates the destroyed handle.
//     /// May be called during serial for_each_handle iteration **ONLY** on the node being iterated.
//     /// May **NOT** be called during parallel for_each_handle iteration.
//     /// May **NOT** be called on the node from which edges are being followed during follow_edges.
//     void destroy_handle(const handle_t& handle);
    
//     /// Create an edge connecting the given handles in the given order and orientations.
//     /// Ignores existing edges.
//     void create_edge(const handle_t& left, const handle_t& right);
    
//     /// Convenient wrapper for create_edge.
//     inline void create_edge(const odgi::edge_t& edge) {
//         create_edge(edge.first, edge.second);
//     }
    
//     /// Remove the edge connecting the given handles in the given order and orientations.
//     /// Ignores nonexistent edges.
//     /// Does not update any stored paths.
//     void destroy_edge(const handle_t& left, const handle_t& right);
    
//     /// Convenient wrapper for destroy_edge.
//     inline void destroy_edge(const edge_t& edge) {
//         destroy_edge(edge.first, edge.second);
//     }
    
//     /// Remove all nodes and edges. Does not update any stored paths.
//     void clear();
    
//     /// Swap the nodes corresponding to the given handles, in the ordering used
//     /// by for_each_handle when looping over the graph. Other handles to the
//     /// nodes being swapped must not be invalidated. If a swap is made while
//     /// for_each_handle is running, it affects the order of the handles
//     /// traversed during the current traversal (so swapping an already seen
//     /// handle to a later handle's position will make the seen handle be visited
//     /// again and the later handle not be visited at all).
//     void swap_handles(const handle_t& a, const handle_t& b);
    
//     /// Alter the node that the given handle corresponds to so the orientation
//     /// indicated by the handle becomes the node's local forward orientation.
//     /// Rewrites all edges pointing to the node and the node's sequence to
//     /// reflect this. Invalidates all handles to the node (including the one
//     /// passed). Returns a new, valid handle to the node in its new forward
//     /// orientation. Note that it is possible for the node's ID to change.
//     /// Does not update any stored paths. May change the ordering of the underlying
//     /// graph.
//     handle_t apply_orientation(const handle_t& handle);
    
//     /// Split a handle's underlying node at the given offsets in the handle's
//     /// orientation. Returns all of the handles to the parts. Other handles to
//     /// the node being split may be invalidated. The split pieces stay in the
//     /// same local forward orientation as the original node, but the returned
//     /// handles come in the order and orientation appropriate for the handle
//     /// passed in.
//     /// Updates stored paths.
//     std::vector<handle_t> divide_handle(const handle_t& handle, const std::vector<size_t>& offsets);
    
//     /// Specialization of divide_handle for a single division point
//     inline std::pair<handle_t, handle_t> divide_handle(const handle_t& handle, size_t offset) {
//         auto parts = divide_handle(handle, std::vector<size_t>{offset});
//         return std::make_pair(parts.front(), parts.back());
//     }

// /**
//  * This is the interface for a handle graph with embedded paths where the paths can be modified.
//  * Note that if the *graph* can also be modified, the implementation will also
//  * need to inherit from MutableHandleGraph, via the combination
//  * MutablePathMutableHandleGraph interface.
//  * TODO: This is a very limited interface at the moment. It will probably need to be extended.
//  */
    
//     /**
//      * Destroy the given path. Invalidates handles to the path and its node steps.
//      */
//     void destroy_path(const path_handle_t& path);

//     /**
//      * Create a path with the given name. The caller must ensure that no path
//      * with the given name exists already, or the behavior is undefined.
//      * Returns a handle to the created empty path. Handles to other paths must
//      * remain valid.
//      */
//     path_handle_t create_path_handle(const std::string& name);
    
//     /**
//      * Append a visit to a node to the given path. Returns a handle to the new
//      * final step on the path which is appended. Handles to prior
//      * steps on the path, and to other paths, must remain valid.
//      */
//     step_handle_t append_step(const path_handle_t& path, const handle_t& to_append);


};
