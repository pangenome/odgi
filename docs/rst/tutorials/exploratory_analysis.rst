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

-------------------------
Build the DRB1-3123 graph
-------------------------

Assuming that your current working directory is the root of the ``odgi`` project, to construct an ``odgi`` file from the
``DRB1-3123`` dataset in ``GFA`` format, execute:

.. code-block:: bash

    odgi build -g test/DRB1-3123.gfa -o DRB1-3123.og

The command creates a file called ``DRB1-3123.og``, which contains the input graph in ``odgi`` format. This graph contains
12 ALT sequences of the `HLA-DRB1 gene <https://www.ncbi.nlm.nih.gov/gene/3123>`_ from the GRCh38 reference genome.

-----------------------------
Visualize the DRB1-3123 graph
-----------------------------

To visualize the graph, execute:

.. code-block:: bash

    odgi viz -i DRB1-3123.og -o DRB1-3123.png -x 500

to obtain the following PNG image:

.. image:: /img/DRB1-3123.png

In this 1-dimensional visualization:

- the graph nodes are arranged from left to right, forming the ``pangenome sequence``;
- the colored bars represent the the paths versus the ``pangenome sequences`` in a binary matrix;
- the path names are visualized on the left;
- the black lines under the paths are the links, which represent the graph topology.

See :ref:`odgi viz` documentation for more information.

----------------------------------
Color respect to the node position
----------------------------------

This is a linearized visualization, but the pangenome graphs are not linear when the embedded genomes present structural
variantion. However, a graph can be optimized for being better visualized in 1-dimension by sorting its nodes properly
(see the :ref:`sorting-layouting` tutorial for more information).

To color the bars respect to the node position in each path, execute:

.. code-block:: bash

    odgi viz -i DRB1-3123.og -o DRB1-3123.du.png -x 500 -du

to obtain the following PNG image:

.. image:: /img/DRB1-3123.du.png

For each path, the brightness goes from light (for the starting position) to black (for the ending position). A linear
genome in a well-sorted graph appears with a smooth brightness gradient in this visualization modality.

Interestingly, the ``>gi|345525392:5000-18402`` path has a brightness gradient which go from right to left. DNA sequence
graphs have two strands, with the node implicitly representing both strands. That gradient indicates that the path is
reversed with respect the ``pangenome sequence``.

------------------------------------
Color respect to the node strandness
------------------------------------

To color the bars respect to the strandness that each node has in each path, execute:

.. code-block:: bash

    odgi viz -i DRB1-3123.og -o DRB1-3123.z.png -x 500 -z

to obtain the following PNG image:

.. image:: /img/DRB1-3123.z.png

A red bar in a path indicates that that region is inverted in that path with respect to the ``pangenome sequence``.


-----------------------------
Build the Lipoprotein A graph
-----------------------------

Assuming that your current working directory is the root of the ``odgi`` project, to construct an ``odgi`` file from the
``LPA`` dataset in ``GFA`` format, execute:

.. code-block:: bash

    odgi build -g test/LPA.gfa -o LPA.og

The command creates a file called ``LPA.og``, which contains the input graph in ``odgi`` format. This graph contains
13 contigs from 7 haploid human genome assemblies from 6 individuals plus the chm13 cell line. The contigs cover the
`Lipoprotein A (LPA) <https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000198670>`_ locus, which encodes the
Apo(a) protein.

-----------------------
Visualize the LPA graph
-----------------------

To visualize the graph, execute:

.. code-block:: bash

    odgi viz -i LPA.og -o LPA.b.png -x 500 -b

to obtain the following PNG image:

.. image:: /img/LPA.b.png

-------------------------------
Color respect to the node depth
-------------------------------

Eukaryotic genomes are characterized by repetitive sequences. These sequences can lead to complex regions in the pangenome
graphs. To identify them, we can analyze the **depth** in the graph. Here we define as **node depth** the number of times
in which the node is crossed by all the paths present in the graph.

To color the bars respect to the mean `depth`, execute:

.. code-block:: bash

    odgi viz -i LPA.og -o LPA.bm.png -x 500 -bm

to obtain the following PNG image:

.. image:: /img/LPA.bm.png

Low depth regions are black, while high depth regions are colored green. Apo(a) proteins vary in size due to a size
polymorphism, the KIV-2 variable numbers of tandem repeats (VNTRs). The VNTR region in the LPA pangenome presents high
`depth`, that becomes evident as a light green stripe in the image.
