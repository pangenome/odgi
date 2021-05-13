.. _sorting-layouting:

#####################
Sorting and Layouting
#####################

========
Synopsis
========

Pangenome graphs built from raw sets of alignments may have complex structures which can introduce difficulty in
downstream analyses, visualization, mapping, and interpretation. Graph sorting aims to find the best node order for
a 1D and 2D layout to simplify these complex regions. Pangenome graphs embed linear pangenomic sequences as paths in
the graph, but to our knowledge, no algorithm takes into account this biological information in the sorting. Moreover,
existing 2D layout methods struggle to deal with large graphs. ``odgi`` implements a new layout algorithm to simplify a pangenome
graph, by using path-guided `stochastic gradient descent <https://ieeexplore.ieee.org/document/8419285>`_
(`PG-SGD <https://docs.google.com/presentation/d/1SfFAtesY6NkSzolo3kN2s3LV5eFunko6KoCv5PkH-YI/edit#slide=id.p>`_) to move a single pair of nodes at a time.
The PG-SGD is memory polite, because it uses a path index, a strict subset of the `xg <https://github.com/vgteam/xg>`_ index. Following a parallelized, lock-free SGD approach,
the PG-SGD can go `Hogwild <https://papers.nips.cc/paper/2011/hash/218a0aefd1d1a4be65601cc6ddc1520e-Abstract.html>`_!
    The 1D path-guided SGD implementation is a key step in general pangenome analyses such as pangenome graph
    linearization and simplification. It is applied in the `PangenomeGraph Builder <https://github.com/pangenome/pggb>`_ (PGGB) pipeline.
This tutorial shows how to sort and visualize a graph in 1D. It explains how to generate a 2D layout of a graph, and how
to take a look at the calculated layout using static and interactive tools.

=====
Steps
=====

----------------------------------
Build the unsorted DRB1-3123 graph
----------------------------------

Assuming that your current working directory is the root of the ``odgi`` project, to construct an ``odgi`` graph from the
``DRB1-3123`` dataset in ``GFA`` format, execute:

.. code-block:: bash

    odgi build -g test/DRB1-3123_unsorted.gfa -o DRB1-3123_unsorted.og

The command creates a file called ``DRB1-3123_unsorted.og``, which contains the input graph in ``odgi`` format. This graph contains
12 ALT sequences of the `HLA-DRB1 gene <https://www.ncbi.nlm.nih.gov/gene/3123>`_ from the GRCh38 reference genome.
The graph was initially created with the `seqwish <https://github.com/ekg/seqwish>`_ variation graph inducer which produces unsorted, raw graphs from
all versus all alignments of the input sequences.

--------------------------------------
Visualize the unsorted DRB1-3123 graph
--------------------------------------

.. code-block:: bash

    odgi viz -i DRB1-3123_unsorted.og -o DRB1-3123_unsorted.png

To obtain the following PNG image:

.. image:: /img/DRB1-3123_unsorted.png

In this 1-Dimensional visualization:

- The graph nodes are arranged from left to right, forming the ``pangenome sequence``.
- The colored bars represent the binned, linearized renderings of the embedded paths versus this ``pangenome sequence`` in a binary matrix.
- The path names are visualized on the left.
- The black lines under the paths are the links, which represent the graph topology.

The graph is very complex and, in this form, the underlying structure is unclear. It is hard to reason from it.
It is a perfect candidate for the 1D PG-SGD: It will find the best 1D node order effectively linearizing it.

--------------------------------------------------------
1D layout metrics of the unsorted DRB1-3123 graph
--------------------------------------------------------

odgi stats

---------------------------------------
Sort the unsorted DRB1-3123 graph in 1D
---------------------------------------

4) odgi sort + a little explanation of what is happening, maybe with a cool GIF too
sort + viz
explain sort update from poster
GIF -> each learning rate update visualized, might help to optimize sorting parameters for given graph
!INFO: not deterministic!

--------------------------------------------------------
2D layout metrics of the unsorted DRB1-3123 graph
--------------------------------------------------------

odgi stats

-----------------------------------------
2D layout of the unsorted DRB1-3123 graph
-----------------------------------------

6) odgi layout + a little explanation of what is happening, maybe with a cool GIF too

--------------------------------------------
Drawing the 2D layout of the DRB1-3123 graph
--------------------------------------------

odgi draw

-----------------------------------------------------------------------------
Interactive visualization of the 2D layout of the DRB1-3123 graph in gfaestus
-----------------------------------------------------------------------------

7) gfaestus
