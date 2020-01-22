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

.. autofunction:: odgi.graph.create_edge

.. autofunction:: odgi.graph.destroy_edge

:class:`odgi.graph`
------------------------------------

.. autofunction:: odgi.graph.apply_ordering

.. autofunction:: odgi.graph.apply_path_ordering

.. autofunction:: odgi.graph.clear

.. autofunction:: odgi.graph.clear_paths

.. autofunction:: odgi.graph.load

.. autofunction:: odgi.graph.optimize

.. autofunction:: odgi.graph.serialize

:class:`odgi.handle`
------------------------------------

.. autofunction:: odgi.graph.apply_orientation

.. autofunction:: odgi.graph.combine_handles

.. autofunction:: odgi.graph.create_handle

.. autofunction:: odgi.graph.destroy_handle

.. autofunction:: odgi.graph.divide_handle

.. autofunction:: odgi.graph.flip

:class:`odgi.path_handle`
------------------------------------

.. autofunction:: odgi.graph.append_step

.. autofunction:: odgi.graph.create_path_handle

.. autofunction:: odgi.graph.destroy_path

.. autofunction:: odgi.graph.prepend_step

.. autofunction:: odgi.graph.set_circularity

:class:`odgi.step_handle`
------------------------------------

.. autofunction:: odgi.graph.rewrite_segment

.. autofunction:: odgi.graph.set_step

.. _accessor:

================
Accessor Methods
================

The following list sorts methods in :class:odgi.graph by what object they return information about.

:class:`odgi.edge`
------------------------------------

.. autofunction:: odgi.graph.edge_handle

:class:`odgi.graph`
------------------------------------

.. autofunction:: odgi.graph.get_node_count

.. autofunction:: odgi.graph.get_path_count

.. autofunction:: odgi.graph.has_node

.. autofunction:: odgi.graph.has_path

.. autofunction:: odgi.graph.max_node_id

.. autofunction:: odgi.graph.min_node_id

.. autofunction:: odgi.graph.to_gfa

:class:`odgi.handle`
------------------------------------

.. autofunction:: odgi.graph.forward

.. autofunction:: odgi.graph.get_degree

.. autofunction:: odgi.graph.get_handle

.. autofunction:: odgi.graph.get_id

.. autofunction:: odgi.graph.get_is_reverse

.. autofunction:: odgi.graph.get_length

.. autofunction:: odgi.graph.get_sequence

.. autofunction:: odgi.graph.get_step_count

.. autofunction:: odgi.graph.has_edge

.. autofunction:: odgi.graph.steps_of_handle
  
:class:`odgi.path_handle`
------------------------------------

.. autofunction:: odgi.graph.get_handle_of_step

.. autofunction:: odgi.graph.get_is_circular

.. autofunction:: odgi.graph.get_path_handle

.. autofunction:: odgi.graph.get_path_name

.. autofunction:: odgi.graph.get_step_count

.. autofunction:: odgi.graph.is_empty

.. autofunction:: odgi.graph.path_back

.. autofunction:: odgi.graph.path_begin

.. autofunction:: odgi.graph.path_end

.. autofunction:: odgi.graph.path_front_end

:class:`odgi.step_handle`
------------------------------------

.. autofunction:: odgi.graph.get_next_step

.. autofunction:: odgi.graph.get_path

.. autofunction:: odgi.graph.get_path_handle_of_step

.. autofunction:: odgi.graph.get_previous_step

.. autofunction:: odgi.graph.has_next_step

.. autofunction:: odgi.graph.has_previous_step

.. autofunction:: odgi.graph.is_path_end

.. autofunction:: odgi.graph.is_path_front_end

.. _iterator:
  
==================
Iteratator Methods
==================

The following list sorts methods in :class:odgi.graph by what kind of iteratee they operate on. 

:class:`odgi.edge`
------------------------------------

.. autofunction:: odgi.graph.follow_edges

:class:`odgi.handle`
------------------------------------

.. autofunction:: odgi.graph.for_each_handle

:class:`odgi.path_handle`
------------------------------------

.. autofunction:: odgi.graph.for_each_path_handle

:class:`odgi.step_handle`
------------------------------------

.. autofunction:: odgi.graph.for_each_step_in_path

.. autofunction:: odgi.graph.for_each_step_on_handle
