.. _glossary:

##########################
Sorted Glossary of Methods
##########################

.. _mutator:

===============
Mutator Methods
===============

The following lists sorts methods in :class:odgi.graph by what types of objects they modify.

:class:`odgi.edge`
------------------------------------

.. function:: create_edge(*args, **kwargs)
   :module: odgi.graph

   Overloaded function.
   
   1. create_edge(self: odgi.graph, arg0: odgi.handle, arg1: odgi.handle) -> None
   
   Create an edge connecting the given handles in the given order and orientations.
   
   2. create_edge(self: odgi.graph, arg0: Tuple[odgi.handle, odgi.handle]) -> None
   
   Create an edge connecting the given handles in the given order and orientations.
   

.. function:: destroy_edge(*args, **kwargs)
   :module: odgi.graph

   Overloaded function.
   
   1. destroy_edge(self: odgi.graph, arg0: Tuple[odgi.handle, odgi.handle]) -> None
   
   Remove the edge connecting the given handles in the given order and orientations.
   
   2. destroy_edge(self: odgi.graph, arg0: odgi.handle, arg1: odgi.handle) -> None
   
   Remove the edge connecting the given handles in the given order and orientations.

:class:`odgi.graph`
------------------------------------

.. function:: apply_ordering(self: odgi.graph, order: List[odgi.handle], compact_ids: bool = False) -> None
   :module: odgi.graph

   Reorder the graph's internal structure to match that given.
   Optionally compact the id space of the graph to match the ordering, from 1->|ordering|.
   

.. function:: apply_path_ordering(self: odgi.graph, arg0: List[odgi.path_handle]) -> None
   :module: odgi.graph

   Reorder the graph's paths as given.
   

.. function:: clear(self: odgi.graph) -> None
   :module: odgi.graph

   Remove all nodes and edges. Does not update any stored paths.
   

.. function:: clear_paths(self: odgi.graph) -> None
   :module: odgi.graph

   Remove all stored paths.
   

.. function:: load(self: odgi.graph, arg0: str) -> None
   :module: odgi.graph

   Load the graph from the given file.
   

.. function:: optimize(self: odgi.graph, allow_id_reassignment: bool = False) -> None
   :module: odgi.graph

   Organize the graph for better performance and memory use.
   

.. function:: serialize(self: odgi.graph, arg0: str) -> None
   :module: odgi.graph

   Save the graph to the given file, returning the number of bytes written.

:class:`odgi.handle`
------------------------------------

.. function:: apply_orientation(self: odgi.graph, arg0: odgi.handle) -> odgi.handle
   :module: odgi.graph

   Alter the node that the given handle corresponds to so the orientation indicated
   by the handle becomes the node's local forward orientation.
   Updates all links and path steps to match the new orientation.
   

.. function:: combine_handles(self: odgi.graph, arg0: List[odgi.handle]) -> odgi.handle
   :module: odgi.graph

   Join handles into a new node, returning the handle of the new node.
   

.. function:: create_handle(*args, **kwargs)
   :module: odgi.graph

   Overloaded function.
   
   1. create_handle(self: odgi.graph, arg0: str) -> odgi.handle
   
   Create a new node with the given sequence and return the handle.
   
   2. create_handle(self: odgi.graph, arg0: str, arg1: int) -> odgi.handle
   
   Create a new node with the given sequence and return the handle.
   

.. function:: destroy_handle(self: odgi.graph, arg0: odgi.handle) -> None
   :module: odgi.graph

   Remove the node belonging to the given handle and all of its edges.
   Does not update any stored paths.
   Invalidates the destroyed handle.
   

.. function:: divide_handle(*args, **kwargs)
   :module: odgi.graph

   Overloaded function.
   
   1. divide_handle(self: odgi.graph, arg0: odgi.handle, arg1: List[int]) -> List[odgi.handle]
   
   Split a handle's underlying node at the given offsets in the handle's orientation.
   Returns the handles to the new parts.
   
   2. divide_handle(self: odgi.graph, arg0: odgi.handle, arg1: int) -> Tuple[odgi.handle, odgi.handle]
   
   Split a handle's underlying node at the given offset in the handle's orientation.
   Returns the handles to the new parts.
   

.. function:: flip(self: odgi.graph, handle: odgi.handle) -> odgi.handle
   :module: odgi.graph

   Flip the handle to the opposite orientation.

:class:`odgi.path_handle`
------------------------------------

.. function:: append_step(self: odgi.graph, arg0: odgi.path_handle, arg1: odgi.handle) -> odgi.step_handle
   :module: odgi.graph

   Append a visit to a node to the given path.
   Returns a handle to the new final step
   on the path which is appended.
   

.. function:: create_path_handle(self: odgi.graph, name: str, is_circular: bool = False) -> odgi.path_handle
   :module: odgi.graph

   Create a path with the given name. The caller must ensure that no path with the
   given name already exists.
   

.. function:: destroy_path(self: odgi.graph, arg0: odgi.path_handle) -> None
   :module: odgi.graph

   Destroy the given path. Invalidates handles to the path and its node steps.
   

.. function:: prepend_step(self: odgi.graph, arg0: odgi.path_handle, arg1: odgi.handle) -> odgi.step_handle
   :module: odgi.graph

   Append a visit to a node to the given path.
   Returns a handle to the new final step on the path which is appended.
   

.. function:: set_circularity(self: odgi.graph, arg0: odgi.path_handle, arg1: bool) -> None
   :module: odgi.graph

   Set if the path is circular or not.

:class:`odgi.step_handle`
------------------------------------

.. function:: rewrite_segment(self: odgi.graph, arg0: odgi.step_handle, arg1: odgi.step_handle, arg2: List[odgi.handle]) -> Tuple[odgi.step_handle, odgi.step_handle]
   :module: odgi.graph

   Replace the path range with the new segment,
   returning the new start and end step handles for the segment.
   

.. function:: set_step(self: odgi.graph, arg0: odgi.step_handle, arg1: odgi.handle) -> odgi.step_handle
   :module: odgi.graph

   Set the step to the given handle, possibly re-linking and cleaning up if needed.

.. _accessor:

================
Accessor Methods
================

The following list sorts methods in :class:odgi.graph by what object they return information about.

:class:`odgi.edge`
------------------------------------

.. function:: edge_handle(self: odgi.graph, arg0: odgi.handle, arg1: odgi.handle) -> Tuple[odgi.handle, odgi.handle]
   :module: odgi.graph

   Return the edge handle for the given pair of handles.

:class:`odgi.graph`
------------------------------------

.. function:: get_node_count(self: odgi.graph) -> int
   :module: odgi.graph

   Return the number of nodes in the graph.
   

.. function:: get_path_count(self: odgi.graph) -> int
   :module: odgi.graph

   Return the path count of the graph
   

.. function:: has_node(self: odgi.graph, node_id: int) -> bool
   :module: odgi.graph

   Return true if the given node is in the graph.
   

.. function:: has_path(self: odgi.graph, arg0: str) -> bool
   :module: odgi.graph

   Return if a path with the givenv name exists in the graph.
   

.. function:: max_node_id(self: odgi.graph) -> int
   :module: odgi.graph

   Return the maximum node id in the graph.
   

.. function:: min_node_id(self: odgi.graph) -> int
   :module: odgi.graph

   Return the minimum node id in the graph.
   

.. function:: to_gfa(self: odgi.graph) -> None
   :module: odgi.graph

   Display as GFA

:class:`odgi.handle`
------------------------------------
  
.. function:: forward(self: odgi.graph, arg0: odgi.handle) -> odgi.handle
   :module: odgi.graph

   Return the forward version of the handle.
   

.. function:: get_degree(self: odgi.graph, arg0: odgi.handle, arg1: bool) -> int
   :module: odgi.graph

   Return the degree of the given node.
   

.. function:: get_handle(self: odgi.graph, node_id: int, is_reverse: bool = False) -> odgi.handle
   :module: odgi.graph

   Return the handle for the given node id.
   

.. function:: get_id(self: odgi.graph, handle: odgi.handle) -> int
   :module: odgi.graph

   Return the id of the given handle.
   

.. function:: get_is_reverse(self: odgi.graph, handle: odgi.handle) -> bool
   :module: odgi.graph

   Return true if the handle refers to the node reverse complement.
   

.. function:: get_length(self: odgi.graph, handle: odgi.handle) -> int
   :module: odgi.graph

   Return the length of the node referred to by the handle.
   

.. function:: get_sequence(self: odgi.graph, handle: odgi.handle) -> str
   :module: odgi.graph


.. function:: get_step_count(*args, **kwargs)
   :module: odgi.graph

   Overloaded function.
   
   1. get_step_count(self: odgi.graph, arg0: odgi.path_handle) -> int
   
   Return the step count of a given path.
   
   2. get_step_count(self: odgi.graph, arg0: odgi.handle) -> int
   
   Return the number of steps on the given handle.
   

.. function:: has_edge(self: odgi.graph, arg0: odgi.handle, arg1: odgi.handle) -> bool
   :module: odgi.graph

   Returns true if the given edge exists
   

.. function:: steps_of_handle(self: odgi.graph, arg0: odgi.handle, arg1: bool) -> List[odgi.step_handle]
   :module: odgi.graph

   Obtain the steps on a given handle.
   

:class:`odgi.path_handle`
------------------------------------

.. function:: get_handle_of_step(self: odgi.graph, arg0: odgi.step_handle) -> odgi.handle
   :module: odgi.graph

   Return the handle that a given step occurs on.
   

.. function:: get_is_circular(self: odgi.graph, arg0: odgi.path_handle) -> bool
   :module: odgi.graph

   Returns true if the path is circular.
   

.. function:: get_path_handle(self: odgi.graph, arg0: str) -> odgi.path_handle
   :module: odgi.graph

   Return the path handle for the named path.
   

.. function:: get_path_name(self: odgi.graph, arg0: odgi.path_handle) -> str
   :module: odgi.graph

   Return the path name for a given path handle.
   

.. function:: get_step_count(*args, **kwargs)
   :module: odgi.graph

   Overloaded function.
   
   1. get_step_count(self: odgi.graph, arg0: odgi.path_handle) -> int
   
   Return the step count of a given path.
   
   2. get_step_count(self: odgi.graph, arg0: odgi.handle) -> int
   
   Return the number of steps on the given handle.
   

.. function:: is_empty(self: odgi.graph, arg0: odgi.path_handle) -> bool
   :module: odgi.graph

   Returns true if the given path is empty, and false otherwise.
   

.. function:: path_back(self: odgi.graph, arg0: odgi.path_handle) -> odgi.step_handle
   :module: odgi.graph

   Return a step handle to the last step, which is arbitrary in the case
   of a circular path.
   

.. function:: path_begin(self: odgi.graph, arg0: odgi.path_handle) -> odgi.step_handle
   :module: odgi.graph

   Return the step handle for the first step in the given path.
   

.. function:: path_end(self: odgi.graph, arg0: odgi.path_handle) -> odgi.step_handle
   :module: odgi.graph

   Return a step handle to a fictitious handle one past the end of the path.
   

.. function:: path_front_end(self: odgi.graph, arg0: odgi.path_handle) -> odgi.step_handle
   :module: odgi.graph

   Return a step handle to a fictitious handle one past the start of the path.

:class:`odgi.step_handle`
------------------------------------

.. function:: get_next_step(self: odgi.graph, arg0: odgi.step_handle) -> odgi.step_handle
   :module: odgi.graph

   Returns a handle to the next step on the path. Calling on an end marker
   step returns the same end marker.
   

.. function:: get_path(self: odgi.graph, arg0: odgi.step_handle) -> odgi.path_handle
   :module: odgi.graph

   Return the path of a given step handle.
   

.. function:: get_path_handle_of_step(self: odgi.graph, arg0: odgi.step_handle) -> odgi.path_handle
   :module: odgi.graph

   Returns a handle to the path that an step is on.
   

.. function:: get_previous_step(self: odgi.graph, arg0: odgi.step_handle) -> odgi.step_handle
   :module: odgi.graph

   Returns a handle to the previous step on the path. Calling on a front
   end marker step returns the same end marker.
   

.. function:: has_next_step(self: odgi.graph, arg0: odgi.step_handle) -> bool
   :module: odgi.graph

   Returns true if the step is not the last step on the path, else false.
   

.. function:: has_previous_step(self: odgi.graph, arg0: odgi.step_handle) -> bool
   :module: odgi.graph

   Returns true if the step is not the first step on the path, else false.
   

.. function:: is_path_end(self: odgi.graph, arg0: odgi.step_handle) -> bool
   :module: odgi.graph

   Returns true if the step handle is an end magic handle.
   

.. function:: is_path_front_end(self: odgi.graph, arg0: odgi.step_handle) -> bool
   :module: odgi.graph

   Returns true if the step handle is a front end magic handle.

.. _iterator:
  
==================
Iteratator Methods
==================

The following list sorts methods in :class:odgi.graph by what kind of iteratee they operate on. 

:class:`odgi.edge`
------------------------------------

.. function:: follow_edges(self: odgi.graph, arg0: odgi.handle, arg1: bool, arg2: Callable[[odgi.handle], bool]) -> bool
   :module: odgi.graph

   Follow edges starting at a given node.

:class:`odgi.handle`
------------------------------------

.. function:: for_each_handle(self: odgi.graph, iteratee: Callable[[odgi.handle], bool], parallel: bool = False) -> bool
   :module: odgi.graph

   Iterate over all the nodes in the graph.

:class:`odgi.path_handle`
------------------------------------

.. function:: for_each_path_handle(self: odgi.graph, arg0: Callable[[odgi.path_handle], bool]) -> bool
   :module: odgi.graph

   Invoke the callback for each path in the graph.

:class:`odgi.step_handle`
------------------------------------

.. function:: for_each_step_in_path(self: odgi.graph, arg0: odgi.path_handle, arg1: Callable[[odgi.step_handle], None]) -> None
   :module: odgi.graph

   Invoke the callback for each step in a given path.
   

.. function:: for_each_step_on_handle(self: odgi.graph, arg0: odgi.handle, arg1: Callable[[odgi.step_handle], bool]) -> bool
   :module: odgi.graph

   Invoke the callback for each of the steps on a given handle.
