.. _api:

##########
python API
##########

..
 automodule:: odgi
   :members:
   commented out.  Use modified autodoc to turn this output into .rst output so readthedocs isn't angry.


.. py:module:: odgi


.. py:class:: edge
   :module: odgi

   edges link two handles together
   

.. py:class:: graph
   :module: odgi

   the odgi graph type
   
   
   .. py:method:: graph.append_step(self: odgi.graph, arg0: odgi.path_handle, arg1: odgi.handle) -> odgi.step_handle
      :module: odgi
   
      Append a visit to a node to the given path.
      Returns a handle to the new final step
      on the path which is appended.
      
   
   .. py:method:: graph.apply_ordering(self: odgi.graph, order: List[odgi.handle], compact_ids: bool = False) -> None
      :module: odgi
   
      Reorder the graph's internal structure to match that given.
      Optionally compact the id space of the graph to match the ordering, from 1->|ordering|.
      
   
   .. py:method:: graph.apply_orientation(self: odgi.graph, arg0: odgi.handle) -> odgi.handle
      :module: odgi
   
      Alter the node that the given handle corresponds to so the orientation indicated
      by the handle becomes the node's local forward orientation.
      Updates all links and path steps to match the new orientation.
      
   
   .. py:method:: graph.apply_path_ordering(self: odgi.graph, arg0: List[odgi.path_handle]) -> None
      :module: odgi
   
      Reorder the graph's paths as given.
      
   
   .. py:method:: graph.clear(self: odgi.graph) -> None
      :module: odgi
   
      Remove all nodes and edges. Does not update any stored paths.
      
   
   .. py:method:: graph.clear_paths(self: odgi.graph) -> None
      :module: odgi
   
      Remove all stored paths.
      
   
   .. py:method:: graph.combine_handles(self: odgi.graph, arg0: List[odgi.handle]) -> odgi.handle
      :module: odgi
   
      Join handles into a new node, returning the handle of the new node.
      
   
   .. py:method:: graph.create_edge(*args, **kwargs)
      :module: odgi
   
      Overloaded function.
      
      1. create_edge(self: odgi.graph, arg0: odgi.handle, arg1: odgi.handle) -> None
      
      Create an edge connecting the given handles in the given order and orientations.
      
      2. create_edge(self: odgi.graph, arg0: Tuple[odgi.handle, odgi.handle]) -> None
      
      Create an edge connecting the given handles in the given order and orientations.
      
   
   .. py:method:: graph.create_handle(*args, **kwargs)
      :module: odgi
   
      Overloaded function.
      
      1. create_handle(self: odgi.graph, arg0: str) -> odgi.handle
      
      Create a new node with the given sequence and return the handle.
      
      2. create_handle(self: odgi.graph, arg0: str, arg1: int) -> odgi.handle
      
      Create a new node with the given sequence and return the handle.
      
   
   .. py:method:: graph.create_path_handle(self: odgi.graph, name: str, is_circular: bool = False) -> odgi.path_handle
      :module: odgi
   
      Create a path with the given name. The caller must ensure that no path with the
      given name already exists.
      
   
   .. py:method:: graph.destroy_edge(*args, **kwargs)
      :module: odgi
   
      Overloaded function.
      
      1. destroy_edge(self: odgi.graph, arg0: Tuple[odgi.handle, odgi.handle]) -> None
      
      Remove the edge connecting the given handles in the given order and orientations.
      
      2. destroy_edge(self: odgi.graph, arg0: odgi.handle, arg1: odgi.handle) -> None
      
      Remove the edge connecting the given handles in the given order and orientations.
      
   
   .. py:method:: graph.destroy_handle(self: odgi.graph, arg0: odgi.handle) -> None
      :module: odgi
   
      Remove the node belonging to the given handle and all of its edges.
      Does not update any stored paths.
      Invalidates the destroyed handle.
      
   
   .. py:method:: graph.destroy_path(self: odgi.graph, arg0: odgi.path_handle) -> None
      :module: odgi
   
      Destroy the given path. Invalidates handles to the path and its node steps.
      
   
   .. py:method:: graph.divide_handle(*args, **kwargs)
      :module: odgi
   
      Overloaded function.
      
      1. divide_handle(self: odgi.graph, arg0: odgi.handle, arg1: List[int]) -> List[odgi.handle]
      
      Split a handle's underlying node at the given offsets in the handle's orientation.
      Returns the handles to the new parts.
      
      2. divide_handle(self: odgi.graph, arg0: odgi.handle, arg1: int) -> Tuple[odgi.handle, odgi.handle]
      
      Split a handle's underlying node at the given offset in the handle's orientation.
      Returns the handles to the new parts.
      
   
   .. py:method:: graph.edge_handle(self: odgi.graph, arg0: odgi.handle, arg1: odgi.handle) -> Tuple[odgi.handle, odgi.handle]
      :module: odgi
   
      Return the edge handle for the given pair of handles.
      
   
   .. py:method:: graph.flip(self: odgi.graph, handle: odgi.handle) -> odgi.handle
      :module: odgi
   
      Flip the handle to the opposite orientation.
      
   
   .. py:method:: graph.follow_edges(self: odgi.graph, arg0: odgi.handle, arg1: bool, arg2: Callable[[odgi.handle], bool]) -> bool
      :module: odgi
   
      Follow edges starting at a given node.
      
   
   .. py:method:: graph.for_each_handle(self: odgi.graph, iteratee: Callable[[odgi.handle], bool], parallel: bool = False) -> bool
      :module: odgi
   
      Iterate over all the nodes in the graph.
      
   
   .. py:method:: graph.for_each_path_handle(self: odgi.graph, arg0: Callable[[odgi.path_handle], bool]) -> bool
      :module: odgi
   
      Invoke the callback for each path in the graph.
      
   
   .. py:method:: graph.for_each_step_in_path(self: odgi.graph, arg0: odgi.path_handle, arg1: Callable[[odgi.step_handle], None]) -> None
      :module: odgi
   
      Invoke the callback for each step in a given path.
      
   
   .. py:method:: graph.for_each_step_on_handle(self: odgi.graph, arg0: odgi.handle, arg1: Callable[[odgi.step_handle], bool]) -> bool
      :module: odgi
   
      Invoke the callback for each of the steps on a given handle.
      
   
   .. py:method:: graph.forward(self: odgi.graph, arg0: odgi.handle) -> odgi.handle
      :module: odgi
   
      Return the forward version of the handle.
      
   
   .. py:method:: graph.get_degree(self: odgi.graph, arg0: odgi.handle, arg1: bool) -> int
      :module: odgi
   
      Return the degree of the given node.
      
   
   .. py:method:: graph.get_handle(self: odgi.graph, node_id: int, is_reverse: bool = False) -> odgi.handle
      :module: odgi
   
      Return the handle for the given node id.
      
   
   .. py:method:: graph.get_handle_of_step(self: odgi.graph, arg0: odgi.step_handle) -> odgi.handle
      :module: odgi
   
      Return the handle that a given step occurs on.
      
   
   .. py:method:: graph.get_id(self: odgi.graph, handle: odgi.handle) -> int
      :module: odgi
   
      Return the id of the given handle.
      
   
   .. py:method:: graph.get_is_circular(self: odgi.graph, arg0: odgi.path_handle) -> bool
      :module: odgi
   
      Returns true if the path is circular.
      
   
   .. py:method:: graph.get_is_reverse(self: odgi.graph, handle: odgi.handle) -> bool
      :module: odgi
   
      Return true if the handle refers to the node reverse complement.
      
   
   .. py:method:: graph.get_length(self: odgi.graph, handle: odgi.handle) -> int
      :module: odgi
   
      Return the length of the node referred to by the handle.
      
   
   .. py:method:: graph.get_next_step(self: odgi.graph, arg0: odgi.step_handle) -> odgi.step_handle
      :module: odgi
   
      Returns a handle to the next step on the path. Calling on an end marker
      step returns the same end marker.
      
   
   .. py:method:: graph.get_node_count(self: odgi.graph) -> int
      :module: odgi
   
      Return the number of nodes in the graph.
      
   
   .. py:method:: graph.get_path(self: odgi.graph, arg0: odgi.step_handle) -> odgi.path_handle
      :module: odgi
   
      Return the path of a given step handle.
      
   
   .. py:method:: graph.get_path_count(self: odgi.graph) -> int
      :module: odgi
   
      Return the path count of the graph
      
   
   .. py:method:: graph.get_path_handle(self: odgi.graph, arg0: str) -> odgi.path_handle
      :module: odgi
   
      Return the path handle for the named path.
      
   
   .. py:method:: graph.get_path_handle_of_step(self: odgi.graph, arg0: odgi.step_handle) -> odgi.path_handle
      :module: odgi
   
      Returns a handle to the path that an step is on.
      
   
   .. py:method:: graph.get_path_name(self: odgi.graph, arg0: odgi.path_handle) -> str
      :module: odgi
   
      Return the path name for a given path handle.
      
   
   .. py:method:: graph.get_previous_step(self: odgi.graph, arg0: odgi.step_handle) -> odgi.step_handle
      :module: odgi
   
      Returns a handle to the previous step on the path. Calling on a front
      end marker step returns the same end marker.
      
   
   .. py:method:: graph.get_sequence(self: odgi.graph, handle: odgi.handle) -> str
      :module: odgi
   
   
   .. py:method:: graph.get_step_count(*args, **kwargs)
      :module: odgi
   
      Overloaded function.
      
      1. get_step_count(self: odgi.graph, arg0: odgi.path_handle) -> int
      
      Return the step count of a given path.
      
      2. get_step_count(self: odgi.graph, arg0: odgi.handle) -> int
      
      Return the number of steps on the given handle.
      
   
   .. py:method:: graph.has_edge(self: odgi.graph, arg0: odgi.handle, arg1: odgi.handle) -> bool
      :module: odgi
   
      Returns true if the given edge exists
      
   
   .. py:method:: graph.has_next_step(self: odgi.graph, arg0: odgi.step_handle) -> bool
      :module: odgi
   
      Returns true if the step is not the last step on the path, else false.
      
   
   .. py:method:: graph.has_node(self: odgi.graph, node_id: int) -> bool
      :module: odgi
   
      Return true if the given node is in the graph.
      
   
   .. py:method:: graph.has_path(self: odgi.graph, arg0: str) -> bool
      :module: odgi
   
      Return if a path with the givenv name exists in the graph.
      
   
   .. py:method:: graph.has_previous_step(self: odgi.graph, arg0: odgi.step_handle) -> bool
      :module: odgi
   
      Returns true if the step is not the first step on the path, else false.
      
   
   .. py:method:: graph.insert_step(self: odgi.graph, arg0: odgi.step_handle, arg1: odgi.step_handle, arg2: odgi.handle) -> odgi.step_handle
      :module: odgi
   
      Insert a visit to a node to the given path between the given steps.
      Returns a handle to the new step on the path which is appended.
      
   
   .. py:method:: graph.is_empty(self: odgi.graph, arg0: odgi.path_handle) -> bool
      :module: odgi
   
      Returns true if the given path is empty, and false otherwise.
      
   
   .. py:method:: graph.is_path_end(self: odgi.graph, arg0: odgi.step_handle) -> bool
      :module: odgi
   
      Returns true if the step handle is an end magic handle.
      
   
   .. py:method:: graph.is_path_front_end(self: odgi.graph, arg0: odgi.step_handle) -> bool
      :module: odgi
   
      Returns true if the step handle is a front end magic handle.
      
   
   .. py:method:: graph.load(self: odgi.graph, arg0: str) -> None
      :module: odgi
   
      Load the graph from the given file.
      
   
   .. py:method:: graph.max_node_id(self: odgi.graph) -> int
      :module: odgi
   
      Return the maximum node id in the graph.
      
   
   .. py:method:: graph.min_node_id(self: odgi.graph) -> int
      :module: odgi
   
      Return the minimum node id in the graph.
      
   
   .. py:method:: graph.optimize(self: odgi.graph, allow_id_reassignment: bool = False) -> None
      :module: odgi
   
      Organize the graph for better performance and memory use.
      
   
   .. py:method:: graph.path_back(self: odgi.graph, arg0: odgi.path_handle) -> odgi.step_handle
      :module: odgi
   
      Return a step handle to the last step, which is arbitrary in the case
      of a circular path.
      
   
   .. py:method:: graph.path_begin(self: odgi.graph, arg0: odgi.path_handle) -> odgi.step_handle
      :module: odgi
   
      Return the step handle for the first step in the given path.
      
   
   .. py:method:: graph.path_end(self: odgi.graph, arg0: odgi.path_handle) -> odgi.step_handle
      :module: odgi
   
      Return a step handle to a fictitious handle one past the end of the path.
      
   
   .. py:method:: graph.path_front_end(self: odgi.graph, arg0: odgi.path_handle) -> odgi.step_handle
      :module: odgi
   
      Return a step handle to a fictitious handle one past the start of the path.
      
   
   .. py:method:: graph.prepend_step(self: odgi.graph, arg0: odgi.path_handle, arg1: odgi.handle) -> odgi.step_handle
      :module: odgi
   
      Append a visit to a node to the given path.
      Returns a handle to the new final step on the path which is appended.
      
   
   .. py:method:: graph.rewrite_segment(self: odgi.graph, arg0: odgi.step_handle, arg1: odgi.step_handle, arg2: List[odgi.handle]) -> Tuple[odgi.step_handle, odgi.step_handle]
      :module: odgi
   
      Replace the path range with the new segment,
      returning the new start and end step handles for the segment.
      
   
   .. py:method:: graph.serialize(self: odgi.graph, arg0: str) -> None
      :module: odgi
   
      Save the graph to the given file, returning the number of bytes written.
      
   
   .. py:method:: graph.set_circularity(self: odgi.graph, arg0: odgi.path_handle, arg1: bool) -> None
      :module: odgi
   
      Set if the path is circular or not.
      
   
   .. py:method:: graph.set_step(self: odgi.graph, arg0: odgi.step_handle, arg1: odgi.handle) -> odgi.step_handle
      :module: odgi
   
      Set the step to the given handle, possibly re-linking and cleaning up if needed.
      
   
   .. py:method:: graph.steps_of_handle(self: odgi.graph, arg0: odgi.handle, arg1: bool) -> List[odgi.step_handle]
      :module: odgi
   
      Obtain the steps on a given handle.
      
   
   .. py:method:: graph.to_gfa(self: odgi.graph) -> None
      :module: odgi
   
      Display as GFA
      

.. py:class:: handle
   :module: odgi

   the handle, which refers to oriented nodes
   

.. py:class:: path_handle
   :module: odgi

   the path handle type, which refers to paths
   

.. py:class:: step_handle
   :module: odgi

   the step handle type, which refers to path paths
