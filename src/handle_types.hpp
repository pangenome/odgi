#ifndef dank_graph_handle_types
#define dank_graph_handle_types
#include <cstdint>
#include <utility>


namespace dsgvg {

/// represents an id
typedef int64_t id_t;
/// represents an offset
typedef std::size_t off_t;
/// represents a position
typedef std::tuple<id_t, bool, off_t> pos_t;
/// represents the internal id of a node traversal
struct handle_t { char data[sizeof(id_t)]; };
/// represents an edge in terms of its endpoints
typedef std::pair<handle_t, handle_t> edge_t;
/// represents the internal id of a path entity
struct path_handle_t { char data[sizeof(int64_t)]; };
/// An occurrence handle is an opaque reference to an occurrence of an oriented node along a path in a graph
/// In dg, it refers to [0], a node id/rank/handle, and [1], a rank within the records on that node
struct occurrence_handle_t { char data[2 * sizeof(int64_t)]; };

}

#endif
