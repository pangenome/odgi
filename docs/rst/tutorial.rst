#########
Tutorial
#########

****************
Creating Graphs
****************
Let's say that you wanted to create the following graph with ODGI:

.. image:: /img/exampleGraph.png

This graph is a combination of nodes (labelled as `n0`, `n1`, ..., `n9`) and directed edges (arrows).

ODGI Objects
=============

Both :class:`odgi.edge` and :class:`odgi.node` are classes that are created and accessed through a :class:`odgi.graph` object.  Individual nodes in the graph are pointed at by :class:`odgi.handle`.

Routes are stored in the graph and accessed through :class:`odgi.path_handle`, which is a series of :class:`odgi.step_handle` linked together.  Each :class:`odgi.step_handle` is points to the node in that step, and also contains

Handles are pointers to specific pieces of the graph, and it is not possible to operate on them directly, aside from comparing whether the objects are equal.  To get information regarding the object that each handle is pointing to, it is necessary to use the corresponding `get` accessor method in :class:`odgi.graph`.

Making a Graph
===============
First, we must create the graph, then make each node and keep track of their handles.

.. code-block:: python

        gr = odgi.graph()
        seq = ["CGA", "TTGG", "CCGT", "C", "GT", "GATAA", "CGG", "ACA", "GCCG", "ATATAAC"]
        n = []
        for s in seq:
                n.append(gr.create_handle(s))

Now we link together these nodes using their handles. Note that each of these handles is directional, so in order to create the bidirectional edge between `n5` and `n8` we use ``create_edge`` twice.

.. code-block:: python

        gr.create_edge(n[0], n[1])
        gr.create_edge(n[1], n[2])
        gr.create_edge(n[2], n[3])
        gr.create_edge(n[2], n[4])
        gr.create_edge(n[3], n[5])
        gr.create_edge(n[5], n[6])
        gr.create_edge(n[5], n[8])
        gr.create_edge(n[6], n[7])
        gr.create_edge(n[6], n[8])
        gr.create_edge(n[7], n[9])
        gr.create_edge(n[8], n[9])
        gr.create_edge(n[8], n[5])


Creating a Path
===============

Say we want to create the following path through this graph:

.. image:: /img/exampleGraphPath.png

To do this, we would need to create a :class:`odgi.path_handle` in the graph, and then append each :class:`odgi.handle` to the end of the path.

.. code-block:: python

        path = gr.create_path_handle("path")
        gr.append_step(path, n[0])
        gr.append_step(path, n[1])
        gr.append_step(path, n[2])
        gr.append_step(path, n[4])
        gr.append_step(path, n[5])
        gr.append_step(path, n[6])
        gr.append_step(path, n[7])
        gr.append_step(path, n[9])

.. warning::

        :func:`odgi.graph.append_step` will not stop you from appending nodes that are not connected to the preceeding node.

.. code-block:: python
        
        # the following code runs without error
        badpath = gr.create_path_handle("badpath")
        gr.append_step(badpath, n[0])
        gr.append_step(badpath, n[3])

Traversing a path
=================

To traverse a path, we need to fetch a series of :class:`odgi.step_handle` from the graph. Note that although we are effectively asking the path for these items in it, all accessor methods are a part of the :class:`odgi.graph` object.

.. code-block:: python

        step = gr.path_begin(path)
        while(gr.has_next_step(step)):
                # get the node handle from the step handle
                current_node_handle = gr.get_handle_of_step(step)
                # ask the node handle for the sequence
                print(gr.get_sequence(current_node_handle))
                # progress to the next step
                step = gr.get_next_step(step)
        current_node_handle = gr.get_handle_of_step(step)
        print(gr.get_sequence(current_node_handle))

Which will output the following:

.. code-block:: 
        
        CGA
        TTGG
        CCGT
        GT
        GATAA
        CGG
        ACA
        ATATAAC

Manipulating a path
===================

.. DANGER::
        Right now none of this works, because insert_step seems to cause a memory leak. 

Say you wanted to edit this path to add the following edges in blue:

.. image:: /img/exampleGraphPath2.png

First, you need to get the step handles corresponding to `n6` and `n7`, and then insert the new nodes to the path with :func:`odgi.graph.insert_step`. *Note that if you had saved the step handles during path creation, it would not be necessary to traverse the path at this step. Decide which objects to save in memory depending on your application*

.. code-block:: python

        step = gr.path_begin(path)
        while(gr.get_handle_of_step(step) != n[6]):
                step = gr.get_next_step(step)
        # now step corresponds to the step handle preceeding our insertion
        next_step = gr.get_next_step(step)
        step = gr.insert_step(step, next_step, n[8]) #1
        step = gr.insert_step(step, next_step, n[5]) #2
        step = gr.insert_step(step, next_step, n[6]) #3

Each call to :func:`odgi.graph.insert_step` returns the step handle pointing to the inserted node.  

During the process of amending the path, the graph looks as follows:

.. figure:: /img/exampleGraphPath.png
        :align: center
        
        Path before additions

.. figure:: /img/exampleGraphPath3.png
        :align: center
        
        Path after line #1

.. figure:: /img/exampleGraphPath4.png
        :align: center
        
        Path after line #2

.. figure:: /img/exampleGraphPath2.png
        :align: center
        
        Path after line #3

*******************************
Saving and Loading ODGI Graphs
*******************************
