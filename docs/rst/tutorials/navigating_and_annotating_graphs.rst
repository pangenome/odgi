################################
Navigating and annotating graphs
################################

========
Synopsis
========

Coordinates in variation graphs are provided by embedded path sequences. The node ids are not meant to be stable.
Translating path positions between graphs reveals the nodes of each graph that are traversed by e.g. the same
reference sequence present in all graphs. This is particularly helpful when two studies make different graphs, but
incorporated the same annotated reference sequence. The graphs may differ, but were well build individually. Picking the
annotated reference sequence as the backbone to compare gene regions, :ref:`odgi position` can help to answer questions
like:

- Where is this gene in the other graph?
- What other paths are crossing it?
- Is there variation going on?

Once a key regions are identified, they can be extracted with :ref:`odgi extract` and visualized in 1D with :ref:`odgi viz` or in
2D with :ref:`odgi layout` and :ref:`odgi draw`.

http://hypervolu.me/~erik/advbioinfo/HPRCy1v2.MHC.fa.gz
http://hypervolu.me/~erik/advbioinfo/HLA_genes.bed
http://hypervolu.me/~erik/advbioinfo/HPRCy1v2.MHC.HLA_genes.bed

=====
Steps
=====

Let's build a small test graph which is used throughout this tutorial:

.. code-block:: bash

    odgi build -g test/k.gfa -o k.og

Using `vg <https://github.com/vgteam/vg>`_, we can obtain the graph in `dot <https://graphviz.org/doc/info/lang.html>`_
format:

.. code-block:: bash

    vg view -F -p -d test/k.gfa > k.gfa.dot

And use Graphviz to visualize the graph:

.. code-block:: bash

    dot -Tpng k.gfa.dot -o k.gfa.dot.png

.. note::
    Ensure that you have a font package for emojis installed. `Noto Color Emoji <https://www.google.com/get/noto/help/emoji/>`_
    is recommended. Or see a list of `Ubuntu packages <https://packages.ubuntu.com/search?keywords=fonts-noto-color-emoji>`_.

.. image:: /img/k.gfa.dot.png

----------------------------------
Path to graph position mapping
----------------------------------

Take a path position in a graph, and display its corresponding graph position.
Let's find out the graph position of the path ``y`` at nucleotide position ``11``.

.. code-block:: bash

	odgi position -i k.og -p y,11,+ -v

Where:

- ``-p`` specifies the path to find the graph position from as a comma-separated triple:

  1. The name of the path.
  2. The nucleotide position of the path.
  3. The orientation at the give nucleotide position of the path.

- ``-v`` ensures that we actually receive graph positions instead of path positions.

We observe on stdout:

.. code-block:: bash

	#source.path.pos	target.graph.pos
	y,11,+	6,1,+

The graph position is encoded as a comma-separated triple: \

1. The node identifier.
2. The nucleotide position of the graph if all nodes where traversed in ascending node identifier order.
3. The orientation of the node.

.. image:: /img/k.gfa.dot_path2graph.png

The red arrow highlights the found graph position.

----------------------------------
Path to path position mapping
----------------------------------

----------------------------------
Graph to path position mapping
----------------------------------

----------------------------------
Offsets in nodes and paths
----------------------------------

----------------------------------
Reference based mapping
----------------------------------

----------------------------------
The lift: Graph to graph position mapping
----------------------------------

Take a path position in a source graph, and use the common paths between the source and target to project it into the target graph.

.. Translate path positions between graphs (odgi position application): we use that to go from a smoothed graph to a
.. consensus graph and vice versa, but we need a more general example of 2 graphs (from different runs, for example).

.. NOTE:
.. - two graphs with different genomes in them except for the reference
.. - two studies make graphs
.. - they are different, but good individually
.. - now let's compare them... I have a variant in some gene I'm interested in in one
.. - where is that in the other graph? what paths are there?
.. - let's pull out the region in both graphs and visualize them


.. 1) Download 2 GFAs from here (?????? and ??????)
.. 2) odgi build + odgi build
.. A) GENERAL EXAMPLE: we need 2 graphs to show a general case (from different runs, for example)
.. 3) ...

.. B) SPECIFIC PGGB EXAMPLE: from consensus graph to smoothed graph
.. odgi is used in productin in pggb (link). Very little explanation, and then explain
.. 4) consensus->smoothed
.. 5 NOT SURE) smoothed->consensus