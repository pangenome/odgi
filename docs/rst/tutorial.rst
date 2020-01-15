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

.. image:: /img/exampleGraphPath.png

*******************************
Saving and Loading ODGI Graphs
*******************************
