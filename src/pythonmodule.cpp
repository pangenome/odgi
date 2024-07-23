// odgi
#include "odgi.hpp"

// Pybind11
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/iostream.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace odgi;

PYBIND11_MODULE(odgi, m)
{

    py::class_<handlegraph::handle_t>(m, "handle", "the handle, which refers to oriented nodes");
    py::class_<handlegraph::path_handle_t>(m, "path_handle", "the path handle type, which refers to paths");

    // see https://github.com/pangenome/odgi/issues/18
    py::class_<handlegraph::step_handle_t>(m, "step_handle", "the step handle type, which refers to path paths")
        .def("path_id", [](handlegraph::step_handle_t &step_handle) {
                return (reinterpret_cast<const int64_t*>(&step_handle)[0]) >> 1;
        })
        .def("is_reverse", [](handlegraph::step_handle_t &step_handle) {
                return (reinterpret_cast<const int64_t*>(&step_handle)[0]) & 1;
        })
        .def("prev_id", [](handlegraph::step_handle_t &step_handle) {
                return (reinterpret_cast<const int64_t*>(&step_handle)[1]);
        })
        .def("prev_rank", [](handlegraph::step_handle_t &step_handle) {
                return (reinterpret_cast<const int64_t*>(&step_handle)[2]);
        })
        .def("next_id", [](handlegraph::step_handle_t &step_handle) {
                return (reinterpret_cast<const int64_t*>(&step_handle)[3]);
        })
        .def("next_rank", [](handlegraph::step_handle_t &step_handle) {
                return (reinterpret_cast<const int64_t*>(&step_handle)[4]);
        })
        .def("__eq__", [](const handlegraph::step_handle_t &self, handlegraph::step_handle_t &step_handle2) {
            for (uint i=0;i<8;i++) {
                if ((reinterpret_cast<const int64_t*>(&self)[i]) != (reinterpret_cast<const int64_t*>(&step_handle2)[i])){
                    return false;
                }
            }
            return true;
        });
    py::class_<handlegraph::edge_t>(m, "edge", "edges link two handles together")
        .def("first", [](const handlegraph::edge_t &edge_handle) {
                return (&edge_handle)->first;
        })
        .def("second", [](const handlegraph::edge_t &edge_handle) {
                return (&edge_handle)->second;
        });
    // Expose class Graph to Python.
    py::class_<odgi::graph_t>(m, "graph", "the odgi graph type")
        .def(py::init())
        .def("has_node",
              &odgi::graph_t::has_node,
             "Return true if the given node is in the graph.",
             py::arg("node_id"))
        .def("get_handle",
             &odgi::graph_t::get_handle,
             "Return the handle for the given node id.",
             py::arg("node_id"),
             py::arg("is_reverse") = false)
        .def("get_id",
             &odgi::graph_t::get_id,
             "Return the id of the given handle.",
             py::arg("handle"))
        .def("get_is_reverse",
             &odgi::graph_t::get_is_reverse,
             "Return true if the handle refers to the node reverse complement.",
             py::arg("handle"))
        .def("flip",
             &odgi::graph_t::flip,
             "Flip the handle to the opposite orientation.",
             py::arg("handle"))
        .def("get_length",
             &odgi::graph_t::get_length,
             "Return the length of the node referred to by the handle.",
             py::arg("handle"))
        .def("get_sequence",
             &odgi::graph_t::get_sequence,
             py::arg("handle"))
         .def("follow_edges",
              [](const odgi::graph_t& g, const handlegraph::handle_t& handle, bool go_left, const std::function<bool(const handlegraph::handle_t&)>& iteratee) {
                  return g.follow_edges(handle, go_left, [&iteratee](const handlegraph::handle_t& h) { iteratee(h); return true; });
              },
             "Loop over all handles to next/previous (False and True, respectively) nodes. Passes them to a callback which returns False to stop iterating and True to continue.  Returns True if we finished and False if we stopped early.")
        .def("for_each_handle",
             [](const odgi::graph_t& g, const std::function<bool(const handlegraph::handle_t&)>& iteratee, bool parallel) {
                 return g.for_each_handle([&iteratee](const handlegraph::handle_t& h){ iteratee(h); return true; }, parallel);
             },
             "Iterate over all the nodes in the graph.",
             py::arg("iteratee"),
             py::arg("parallel") = false)
        .def("get_node_count",
             &odgi::graph_t::get_node_count,
             "Return the number of nodes in the graph.")
        .def("min_node_id",
             &odgi::graph_t::min_node_id,
             "Return the minimum node id in the graph.")
        .def("max_node_id",
             &odgi::graph_t::max_node_id,
             "Return the maximum node id in the graph.")
        .def("set_id_increment ",
             &odgi::graph_t::set_id_increment,
             "Set a base increment for the node id space.")
        .def("get_degree",
             &odgi::graph_t::get_degree,
             "Return the degree of the given node.")
        .def("forward",
             &odgi::graph_t::forward,
             "Return the forward version of the handle.")
        .def("edge_handle",
             &odgi::graph_t::edge_handle,
             "Return the edge handle for the given pair of handles.")
        .def("has_path",
             &odgi::graph_t::has_path,
             "Return if a path with the given name exists in the graph.")
        .def("get_path_handle",
             &odgi::graph_t::get_path_handle,
             "Return the path handle for the named path.")
        .def("get_path_name",
             &odgi::graph_t::get_path_name,
             "Return the path name for a given path handle.")
        .def("get_step_count",
             [](const odgi::graph_t& g, const handlegraph::path_handle_t& path_handle) { return g.get_step_count(path_handle); },
             "Return the step count of a given path.")
        .def("get_path_count",
             &odgi::graph_t::get_path_count,
             "Return the path count of the graph")
        .def("steps_of_handle",
             &odgi::graph_t::steps_of_handle,
             "Obtain the steps on a given handle.")
        .def("for_each_path_handle",
             [](const odgi::graph_t& g,
                const std::function<bool(const handlegraph::path_handle_t&)>& iteratee) {
                 return g.for_each_path_handle([&iteratee](const handlegraph::path_handle_t& p) {
                         iteratee(p); return true;
                     });
             },
             "Invoke the callback for each path in the graph.")
        .def("for_each_step_on_handle",
             [](const odgi::graph_t& g,
                const handlegraph::handle_t& handle,
                const std::function<bool(const handlegraph::step_handle_t&)>& iteratee) {
                 return g.for_each_step_on_handle(handle, [&iteratee](const handlegraph::step_handle_t& s) {
                         iteratee(s); return true;
                     });
             },
             "Invoke the callback for each of the steps on a given handle.")
        .def("get_step_count",
             [](const odgi::graph_t& g, const handlegraph::handle_t& handle) { return g.get_step_count(handle); },
             "Return the number of steps on the given handle.")
        .def("get_handle_of_step",
             &odgi::graph_t::get_handle_of_step,
             "Return the handle that a given step occurs on.")
        .def("get_path",
             &odgi::graph_t::get_path,
             "Return the path of a given step handle.")
        .def("path_begin",
             &odgi::graph_t::path_begin,
             "Return the step handle for the first step in the given path.")
        .def("path_end",
             &odgi::graph_t::path_end,
             "Return a step handle to a fictitious handle one past the end of the path.")
        .def("path_back",
             &odgi::graph_t::path_back,
             "Return a step handle to the last step, which is arbitrary in the case\nof a circular path.")
        .def("path_front_end",
             &odgi::graph_t::path_front_end,
             "Return a step handle to a fictitious handle one past the start of the path.")
        .def("is_path_front_end",
             &odgi::graph_t::is_path_front_end,
             "Returns true if the step handle is a front end magic handle.")
        .def("is_path_end",
             &odgi::graph_t::is_path_end,
             "Returns true if the step handle is an end magic handle.")
        .def("has_next_step",
             &odgi::graph_t::has_next_step,
             "Returns true if the step is not the last step on the path, else false.")
        .def("has_previous_step",
             &odgi::graph_t::has_previous_step,
             "Returns true if the step is not the first step on the path, else false.")
        .def("get_next_step",
             &odgi::graph_t::get_next_step,
             "Returns a handle to the next step on the path. Calling on an end marker\nstep returns the same end marker.")
        .def("get_previous_step",
             &odgi::graph_t::get_previous_step,
             "Returns a handle to the previous step on the path. Calling on a front\nend marker step returns the same end marker.")
        .def("get_path_handle_of_step",
             &graph_t::get_path_handle_of_step,
             "Returns a handle to the path that an step is on.")
        /*
        .def("get_ordinal_rank_of_step",
             &odgi::graph_t::get_ordinal_rank_of_step,
             "Returns the 0-based ordinal rank of a step on a path. (warning: not implemented in odgi)")
        */
        .def("is_empty",
             &odgi::graph_t::is_empty,
             "Returns true if the given path is empty, and false otherwise.")
        .def("for_each_step_in_path",
             &odgi::graph_t::for_each_step_in_path,
             "Invoke the callback for each step in a given path.")
        .def("get_is_circular",
             &odgi::graph_t::get_is_circular,
             "Returns true if the path is circular.")
        .def("set_circularity",
             &odgi::graph_t::set_circularity,
             "Set if the path is circular or not.")
        .def("create_handle",
             [](odgi::graph_t& g, const std::string& sequence) { return g.create_handle(sequence); },
             "Create a new node with the given sequence and return the handle.")
        .def("create_handle",
             [](odgi::graph_t& g, const std::string& sequence, const handlegraph::nid_t& id) { return g.create_handle(sequence, id); },
             "Create a new node with the given sequence and return the handle.")
        .def("destroy_handle",
             &odgi::graph_t::destroy_handle,
             "Remove the node belonging to the given handle and all of its edges.\nDoes not update any stored paths.\nInvalidates the destroyed handle.")
        .def("create_edge",
             [](odgi::graph_t& g, const handlegraph::handle_t& from, const handlegraph::handle_t& to) { return g.create_edge(from, to); },
             "Create an edge connecting the given handles in the given order and orientations.")
        .def("create_edge",
             [](odgi::graph_t& g, const handlegraph::edge_t& edge) { return g.create_edge(edge); },
             "Create an edge connecting the given handles in the given order and orientations.")
        .def("has_edge",
             &odgi::graph_t::has_edge,
             "Returns true if the given edge exists")
        .def("destroy_edge",
             [](odgi::graph_t& g, const handlegraph::edge_t& edge) { return g.destroy_edge(edge); },
             "Remove the edge connecting the given handles in the given order and orientations.")
        .def("destroy_edge",
             [](odgi::graph_t& g, const handlegraph::handle_t& from, const handlegraph::handle_t& to) { return g.destroy_edge(from, to); },
             "Remove the edge connecting the given handles in the given order and orientations.")
        .def("clear",
             &odgi::graph_t::clear,
             "Remove all nodes and edges. Does not update any stored paths.")
        .def("clear_paths",
             &odgi::graph_t::clear_paths,
             "Remove all stored paths.")
        .def("apply_ordering",
             &odgi::graph_t::apply_ordering,
             "Reorder the graph's internal structure to match that given.\nOptionally compact the id space of the graph to match the ordering, from 1->|ordering|.",
             py::arg("order"),
             py::arg("compact_ids") = false)
        .def("optimize",
             &odgi::graph_t::optimize,
             "Organize the graph for better performance and memory use.",
             py::arg("allow_id_reassignment") = false)
        .def("apply_path_ordering",
             &odgi::graph_t::apply_path_ordering,
             "Reorder the graph's paths as given.")
        .def("apply_orientation",
             &odgi::graph_t::apply_orientation,
             "Alter the node that the given handle corresponds to so the orientation indicated\nby the handle becomes the node's local forward orientation.\nUpdates all links and path steps to match the new orientation.")
        .def("divide_handle",
             [](odgi::graph_t& g, const handlegraph::handle_t& handle, const std::vector<size_t>& offsets) {
                 return g.divide_handle(handle, offsets);
             },
             "Split a handle's underlying node at the given offsets in the handle's orientation.\nReturns the handles to the new parts.")
        .def("divide_handle",
             [](odgi::graph_t& g, const handlegraph::handle_t& handle, size_t offset) {
                 return g.divide_handle(handle, offset);
             },
             "Split a handle's underlying node at the given offset in the handle's orientation.\nReturns the handles to the new parts.")
        .def("combine_handles",
             &odgi::graph_t::combine_handles,
             "Join handles into a new node, returning the handle of the new node.")
        .def("destroy_path",
             &odgi::graph_t::destroy_path,
             "Destroy the given path. Invalidates handles to the path and its node steps.")
        .def("create_path_handle",
             &odgi::graph_t::create_path_handle,
             "Create a path with the given name. The caller must ensure that no path with the\ngiven name already exists.",
             py::arg("name"),
             py::arg("is_circular") = false)
        .def("prepend_step",
             &odgi::graph_t::prepend_step,
             "Append a visit to a node to the given path.\nReturns a handle to the new final step on the path which is appended.")
        .def("append_step",
             &odgi::graph_t::append_step,
             "Append a visit to a node to the given path.\nReturns a handle to the new final step\non the path which is appended.")
        .def("insert_step",
             &odgi::graph_t::insert_step,
             "Insert a visit to a node to the given path between the given steps.\nReturns a handle to the new step on the path which is appended.")
        .def("set_step",
             &odgi::graph_t::set_step,
             "Set the step to the given handle, possibly re-linking and cleaning up if needed.")
        .def("rewrite_segment",
             &odgi::graph_t::rewrite_segment,
             "Replace the path range with the new segment,\nreturning the new start and end step handles for the segment.")
        /*
        .def("display",
             &odgi::graph_t::display,
             "A helper function to visualize the state of the graph")
        */
        .def("to_gfa",
             [](const odgi::graph_t& g) {
                 py::scoped_ostream_redirect stream(
                     std::cout,
                     py::module::import("sys").attr("stdout")
                     );
                 g.to_gfa(std::cout);
             },
             "Display as GFA")
        .def("serialize",
             [](odgi::graph_t& g, const std::string& file) {
                 std::ofstream out(file.c_str());
                 g.serialize(out);
             },
             "Save the graph to the given file, returning the number of bytes written.")
        .def("load",
             [](odgi::graph_t& g, const std::string& file) {
                 std::ifstream in(file.c_str());
                 g.deserialize(in);
             },
             "Load the graph from the given file.")
        // Definition of class_<odgi::graph_t> ends here.
    ;

}
