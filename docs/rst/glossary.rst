.. _glossary:

##########################
Sorted Glossary of Methods
##########################

.. _mutator:

===============
Mutator Methods
===============

The following lists sorts methods in :class:`odgi.graph` by what types of objects they modify.

:class:`odgi.edge`

* :func:`odgi.graph.create_edge`

* :func:`odgi.graph.destroy_edge`

:class:`odgi.graph`

* :func:`odgi.graph.apply_ordering`

* :func:`odgi.graph.apply_path_ordering`

* :func:`odgi.graph.clear`

* :func:`odgi.graph.clear_paths`

* :func:`odgi.graph.load`

* :func:`odgi.graph.optimize`

* :func:`odgi.graph.serialize`

:class:`odgi.handle`

* :func:`odgi.graph.apply_orientation`

* :func:`odgi.graph.combine_handles`

* :func:`odgi.graph.create_handle`

* :func:`odgi.graph.destroy_handle`

* :func:`odgi.graph.divide_handle`

* :func:`odgi.graph.flip`

:class:`odgi.path_handle`

* :func:`odgi.graph.append_step`

* :func:`odgi.graph.create_path_handle`

* :func:`odgi.graph.destroy_path`

* :func:`odgi.graph.prepend_step`

* :func:`odgi.graph.set_circularity`

:class:`odgi.step_handle`

* :func:`odgi.graph.rewrite_segment`

* :func:`odgi.graph.set_step`

.. _accessor:

================
Accessor Methods
================

The following list sorts methods in :class:`odgi.graph` by what object they return information about.

:class:`odgi.edge`

* :func:`odgi.graph.edge_handle`

:class:`odgi.graph`

* :func:`odgi.graph.get_node_count`

* :func:`odgi.graph.get_path_count`

* :func:`odgi.graph.has_node`

* :func:`odgi.graph.has_path`

* :func:`odgi.graph.max_node_id`

* :func:`odgi.graph.min_node_id`

* :func:`odgi.graph.to_gfa`

:class:`odgi.handle`

* :func:`odgi.graph.forward`

* :func:`odgi.graph.get_degree`

* :func:`odgi.graph.get_handle`

* :func:`odgi.graph.get_id`

* :func:`odgi.graph.get_is_reverse`

* :func:`odgi.graph.get_length`

* :func:`odgi.graph.get_sequence`

* :func:`odgi.graph.get_step_count`

* :func:`odgi.graph.has_edge`

* :func:`odgi.graph.steps_of_handle`
  
:class:`odgi.path_handle`

* :func:`odgi.graph.get_handle_of_step`

* :func:`odgi.graph.get_is_circular`

* :func:`odgi.graph.get_path_handle`

* :func:`odgi.graph.get_path_name`

* :func:`odgi.graph.get_step_count`

* :func:`odgi.graph.is_empty`

* :func:`odgi.graph.path_back`

* :func:`odgi.graph.path_begin`

* :func:`odgi.graph.path_end`

* :func:`odgi.graph.path_front_end`

:class:`odgi.step_handle`

* :func:`odgi.graph.get_next_step`

* :func:`odgi.graph.get_path`

* :func:`odgi.graph.get_path_handle_of_step`

* :func:`odgi.graph.get_previous_step`

* :func:`odgi.graph.has_next_step`

* :func:`odgi.graph.has_previous_step`

* :func:`odgi.graph.is_path_end`

* :func:`odgi.graph.is_path_front_end`

.. _iterator:
  
==================
Iteratator Methods
==================

The following list sorts methods in :class:`odgi.graph` by what kind of iteratee they operate on. 

:class:`odgi.edge`

* :func:`odgi.graph.follow_edges`

:class:`odgi.graph`

:class:`odgi.handle`

* :func:`odgi.graph.for_each_handle`

:class:`odgi.path_handle`

* :func:`odgi.graph.for_each_path_handle`

:class:`odgi.step_handle`

* :func:`odgi.graph.for_each_step_in_path`

* :func:`odgi.graph.for_each_step_on_handle`
