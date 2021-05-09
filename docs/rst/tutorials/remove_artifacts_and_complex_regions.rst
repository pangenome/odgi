####################################
Remove artifacts and complex regions
####################################

========
Synopsis
========

De novo assemblies may present errors (mis-assembled sequences, e.g., misjoins and erroneous insertions/deletions) or
sequences that are difficult to align (e.g., centromeres). These issues lead to pangenome graphs with artifacts and/or
very complex regions. To make the downstream analyses easier, as read mapping against the graph or graph visualization,
pangenome graphs can be simplified by applying a set of ``odgi`` tools.

.. note::
   This is an advanced tutorial. It is recommended that you follow the other tutorials before tackling this one.


========
XXXXXXXX
========

.. code-block:: bash

.. TODO: do not put a graph with the consensus paths, to simplify the tutorial

Download the pangenome graph of the `Human chromosome 8 <xxx>`_ in ``GFA`` format, and convert it to a graph in ``odgi``
format:

.. code-block:: bash

    odgi build -g chr8.pan.gfa -o chr8.pan.og

The command creates a file called ``chr8.pan.og``, which contains the input graph in ``odgi`` format.

Identify the regions in the pangenome graph with low coverage, which can present artifacts:

.. code-block:: bash

    odgi depth -i chr8.pan.og -w 100:0:1 > chr8.pan.low_coverage.bed

The ``chr8.pan.low_coverage.bed`` file reports the regions with coverage between 0 and 1. Regions closer to 100 bp have
been merged into a single region. The file is in ``BED`` format.

Identify the regions in the pangenome graph with high coverage, which indicate complex regions:

.. code-block:: bash

    odgi depth -i chr8.pan.og -w 100000:5000:100000000 > chr8.pan.high_coverage.bed

The ``chr8.pan.high_coverage.bed`` file reports the regions with coverage between 5000 and 100000000. Regions closer to
100000 bp have been merged into a single region. The file is in ``BED`` format.

Put all regions larger than 10000 bps in the same ``BED`` format file, merging the adjacent ranges:

.. code-block:: bash

    (awk '$3 - $2 > '10000 chr8.pan.low_coverage.bed ; \
        awk '$3 - $2 > '100000 chr8.pan.high_coverage.bed ) | \
        bedtools sort | bedtools merge > regions_to_remove.bed

Little regions are discarded to avoid too fragmented cleaned graphs.

Clean the pangenome graphs, removing the identified regions:

.. code-block:: bash

    odgi paths -i chr8.pan.og -L | grep '^grch38#chr8' > path_to_fully_retain.txt

    odgi extract -i chr8.pan.og -P \
         --inverse \
         -b regions_to_remove.bed \
         -R path_to_fully_retain.txt -o - | \
         odgi explode -i - -p chr8.pan.clean.exp -b 1 -s P -O

The ``path_to_fully_retain.txt`` contains the name of the reference to fully retain in the resulting graph. The
``--inverse`` specifies to remove the regions specifyed in the ``regions_to_remove.bed`` file.

To have basic information on the cleaned graph, execute:

.. code-block:: bash

    odgi stats -i chr8.pan.clean.exp.8.og -S | column -t

.. code-block:: none

    #length    nodes    edges    paths
    149046153  4044095  5600776  65354

To visualize the cleaned graph, first sort it:

.. code-block:: bash

    odgi sort -p Y -i chr8.pan.clean.exp.8.og -o chr8.pan.clean.og -P

and then execute:

.. code-block:: bash

    odgi paths -i chr8.pan.og -L | cut -f 1,2 -d '#' | uniq > haplotype_names.txt

    odgi viz -i chr8.pan.clean.og -x 1000 -o chr8.pan.clean.png -M haplotype_names.txt

to obtain the following PNG image:

.. image:: /img/chr8.pan.clean.png

This 1-dimensional visualization shows that all centromeres have been removed. Indeed, they present high coverage being
very complex regions. Only the GRCh38 reference centromere is present because it was explicitly preserved during the
removal phase of the low and high coverage regions.

Moreover, for two haplotypes (xxx and xxx), a region close to their
centromere is erroneously absent. This may be due to an under-alignment, which led to the generation of low-coverage nodes
in the pangenome graph, which were removed during the removal phase.
