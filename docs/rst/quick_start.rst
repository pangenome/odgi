===========
Quick Start
===========

.. toctree::
    :maxdepth: 1

A pangenome models the full set of genomic elements in a given species or clade. It can efficiently be encoded in the
form of a variation graph, a type of sequence graph that embeds the linear sequences as paths in the graphs themselves.

To exchange pangenomes, the community frequently uses a subset of version 1 of the ``GFA`` (Graphical Fragment Assembly)
format (`GFAv1 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8006571/#FN8>`_). However, ``odgi`` works on a dynamic
succinct variation graph representation, the ``odgi`` format.

Assuming that your current working directory is the root of the ``odgi`` project, to construct an ``odgi`` file from a
``GFA`` file, execute:

.. code-block:: bash

    odgi build -g test/DRB1-3123.gfa -o DRB1-3123.og

The command creates a file called ``DRB1-3123.og``, which contains the input graph in ``odgi`` format.

To have basic information on the graph, execute:

.. code-block:: bash

    odgi stats -i DRB1-3123.og -S

.. code-block:: none

    #length nodes   edges   paths
    21997   4955    6777    12

This graph file has the following properties:

    - the sum of the lengths of all its nodes is equal to 21997 nucleotides;
    - it has 4955 nodes, 6777 edges, and 12 paths.

The path's names are:

.. code-block:: bash

    odgi paths -i DRB1-3123.og -L

.. code-block:: none

    gi|568815592:32578768-32589835
    gi|568815529:3998044-4011446
    gi|568815551:3814534-3830133
    gi|568815561:3988942-4004531
    gi|568815567:3779003-3792415
    gi|568815569:3979127-3993865
    gi|345525392:5000-18402
    gi|29124352:124254-137656
    gi|28212469:126036-137103
    gi|28212470:131613-146345
    gi|528476637:32549024-32560088
    gi|157702218:147985-163915

We can obtain their sequences in ``FASTA`` format:

.. code-block:: bash

    odgi paths -i DRB1-3123.og -f > paths.fasta
    head paths.fasta -n 2

.. code-block:: none

    >gi|568815592:32578768-32589835
    ATTTAACTCCATCTTTGAGAAACATTTAATAATGTAATGTGTTTGTCATACAGGGTGAATACAGATGCACGGG...

To visualize the graph, execute

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


