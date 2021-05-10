####################
Exploratory analysis
####################

========
Synopsis
========

Visualizing the pangenome graphs we can gain insight into the mutual relationship between the embedded genomes and their
variation. However, complex, nonlinear graph structures are difficult to present in a convenient number of dimensions.
``odgi`` offers a binned, linearized 1-dimensional rendering able to display gigabase scale pangenomes, and several
visualization modalities to grasp the complexity in the graphs.


=====
Steps
=====

.. TODO: to find a dataset

Download the pangenome graph of the `xxxx <xxx>`_ in ``GFA`` format, and convert it to a graph in ``odgi``
format:

.. code-block:: bash

    odgi build -g xxx.gfa -o xxx.og

The command creates a file called ``xxx.og``, which contains the input graph in ``odgi`` format.

To visualize the graph, execute:

.. code-block:: bash

    odgi viz -i xxx.og -o xxx.png -x 1000

to obtain the following PNG image:

.. image:: /img/xxx.png

In this 1-dimensional visualization:

- the graph nodes are arranged from left to right, forming the ``pangenome sequence``;
- the colored bars represent the the paths versus the ``pangenome sequences`` in a binary matrix;
- the path names are visualized on the left;
- the black lines under the paths are the links, which represent the graph topology.

This is a linearized visualization, but the pangenome graphs are not linear when the embedded genomes present structural
variantion. However, a graph can be optimized for being better visualized in 1-dimension by sorting its nodes properly
(see :ref:`XXXsorting-stuffXXX` documentation for more information). To color the bars respect to the node position in
each path, execute:

.. code-block:: bash

    odgi viz -i xxx.og -o xxx.png -x 1000 -du

to obtain the following PNG image:

.. image:: /img/xxx.png

For each path, the brightness goes from light (for the starting position) to black (for the ending position). A linear
genome in a well-sorted graph appears with a smooth brightness gradient in this visualization modality.

DNA sequence graphs have two strands, with the node implicitly representing both strands. To color the bars respect to
the strandness that each node has in each path, execute:

.. code-block:: bash

    odgi viz -i xxx.og -o xxx.png -x 1000 -z

to obtain the following PNG image:

.. image:: /img/xxx.png

A red bar in a path indicates that that region is inverted in that path with respect to the ``pangenome sequence``.

Eukaryotic genomes are characterized by repetitive sequences. These sequences can lead to complex regions in the pangenome
graphs. To identify them, we can analyze the **depth** in the graph. Here we define as **node depth** the number of times
in which the node is crossed by all the paths present in the graph.

To color the bars respect to the mean ``depth`, execute:

.. code-block:: bash

    odgi viz -i xxx.og -o xxx.png -x 1000 -bm

to obtain the following PNG image:

.. image:: /img/xxx.png

The repetitive regions in the XXX pangenome present high ``depth`, that is evident as a light green stripe in the image.
