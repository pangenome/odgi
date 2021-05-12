######################
Extract selected `loci`
######################

========
Synopsis
========

Pangenome graphs model the full set of genomic elements in a gives species or clade. Nevertheless, downstream analyses
may require focusing on specific elements. ``odgi`` allows to extract `loci` of interest from the pangenome graph, for
easily working on them.


=====
Steps
=====


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

.. -----------------------------
..  Visualize the DRB1-3123 graph
..  -----------------------------

.. To visualize the graph, execute:

.. .. code-block:: bash

..     odgi layout -i DRB1-3123.og -o DRB1-3123.lay -P

..     odgi draw -i DRB1-3123.og -c DRB1-3123.lay -p DRB1-3123.png


.. to obtain the following PNG image:

.. .. image:: /img/DRB1-3123.draw.png

.. This 2-dimensional visualization shows the graph topology, where each black line representing a node.

------------------
Display path names
------------------

The path's names are:

.. code-block:: bash

    odgi paths -i LPA.og -L

.. code-block:: none

    chm13__LPA__tig00000001
    HG002__LPA__tig00000001
    HG002__LPA__tig00000005
    HG00733__LPA__tig00000001
    HG00733__LPA__tig00000008
    HG01358__LPA__tig00000002
    HG01358__LPA__tig00000010
    HG02572__LPA__tig00000005
    HG02572__LPA__tig00000001
    NA19239__LPA__tig00000002
    NA19239__LPA__tig00000006
    NA19240__LPA__tig00000001
    NA19240__LPA__tig00000012

---------------------------------
Obtain the Lipoprotein A variants
---------------------------------

To see variants for the two contigs of the ``HG02572`` sample, execute:

.. code-block:: bash

     zgrep '^##' test/LPA.chm13__LPA__tig00000001.vcf.gz -v | head -n 5 | cut -f 1-9,16,17 | column -t

.. code-block:: none

    #CHROM                   POS  ID      REF  ALT  QUAL  FILTER  INFO                        FORMAT  HG02572__LPA__tig00000001  HG02572__LPA__tig00000005
    chm13__LPA__tig00000001  6    >3>6    T    A    60    .       AC=1;AF=0.5;AN=2;LV=0;NS=2  GT      1                          0
    chm13__LPA__tig00000001  134  >6>8    G    GT   60    .       AC=1;AF=0.5;AN=2;LV=0;NS=2  GT      1                          0
    chm13__LPA__tig00000001  156  >8>11   A    G    60    .       AC=1;AF=0.5;AN=2;LV=0;NS=2  GT      0                          1
    chm13__LPA__tig00000001  252  >11>13  GT   G    60    .       AC=2;AF=1;AN=2;LV=0;NS=2    GT      1                          1

The pangenome graphs embed all the mutual relationship between the embedded genomes and their variation. In this example,
the variants are called respect to the ``chm13__LPA__tig00000001`` contig, which was used as reference path.

The insertion at position 136 (G > GT) is present in only one of the  ``HG02572``'s contig (``HG02572__LPA__tig00000001``).
To extract the sub-graph where this insertion falls, execute:

.. code-block:: bash
    odgi extract -i LPA.og -n 23 -c 1 -o LPA.21_23_G_GT.og

The instruction extracts:
- the node with ID 23 (``-n 23``),
- the nodes reachable from this node following a single edge (`-c 1`) in the graph topology,
- the edges connecting all the extracted nodes, and
- the paths traversing all the extracted nodes.

To have basic information on the sub-graph, execute:

.. code-block:: bash

    odgi stats -i LPA.21_23_G_GT.og -S

.. code-block:: none

    #length nodes   edges   paths
    644     5       6       3

The extracted path's names are:

.. code-block:: bash

    odgi paths -i LPA.21_23_G_GT.og -L

.. code-block:: none

    chm13__LPA__tig00000001:997-1640
    HG02572__LPA__tig00000005:999-1641
    HG02572__LPA__tig00000001:1035-1678

The sub-graph contains the contig used as reference in the ``VCF`` file, and the two ``HG02572``'s contigs.

-----------------------
Visualize the sub-graph
-----------------------

To visualize the sub-graph, we can also use external tools as `Bandage <https://github.com/rrwick/Bandage>`_, which
supports grpahs in ``GFA`` format. To covert the graph in ``odgi`` format in a graph in ``GFA`` format, execute:

.. code-block:: bash
    odgi view -i LPA.21_23_G_GT.og -g > LPA.21_23_G_GT.gfa

Then, open the ``LPA.21_23_G_GT.gfa`` file with ``Bandage``.

.. image:: /img/LPA.21_23_G_GT.png

The image shows the graph topology, where each colored rectangle represents a node. In particular, three paths support
nodes with ID 21 and 23, and only one path supports the node with ID 22. The node with ID 22 represents in the graph the
additional nucleotide ``T`` presents in the ``HG02572__LPA__tig00000001`` contig as an insertion.

--------------------------
Get the Human chr8 dataset
--------------------------

Download the pangenome graph of the `Human chromosome 8 <https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2021_05_06_pggb/gfas/chr8.pan.gfa.gz>`_
in ``GFA`` format, decompress it, and convert it to a graph in ``odgi`` format:

.. code-block:: bash

    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2021_05_06_pggb/gfas/chr8.pan.gfa.gz
    gunzip chr8.pan.gfa.gz

    odgi build -g chr8.pan.gfa -o chr8.pan.og

The last command creates a file called ``chr8.pan.og``, which contains the input graph in ``odgi`` format. This graph contains
88 haploid, phased human genome assemblies from 44 individuals, plus the chm13 and GRCh38 reference genomes.

-----------------------
Extract the centromeres/beta-defensin cluster
-----------------------


-----------------------
DO SOMETHING
-----------------------

8) odgi viz/layout/bandage
