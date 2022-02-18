.. _faqs:

####
FAQs
####

.. toctree::
    :maxdepth: 1

How can I import reads from a FASTQ or FASTA file into an existing graph?
=========================================================================

``odgi`` can work with any kind of graphs in the GFAv1 format,
including graphs constructed with e.g. ``vg construct`` containing short read data. However, ``odgi`` does not construct nor
extend existing graphs. But, there already are other specialised tools to integrate such sequences:

1. Graph constructed from long read or sequence data, extension with short reads or sequences
---------------------------------------------------------------------------------------------

As the graph might be very complex, a necessary step might be to prune the graph with ``vg prune``. Or removing complex
regions following the ODGI tutorial :ref:`remove_artifacts_and_complex_regions`. These steps might be necessary in order
to build the indices required for `Giraffe <https://www.science.org/doi/epdf/10.1126/science.abg8871>`_. Then map the
sequences to the graph with ``vg giraffe``. The resulting `GAM <https://github.com/vgteam/vg/wiki/File-Formats#gam-graph-alignment--map-vgs-bam>`_ file
can be used with ``vg augment`` to extend the existing graph with the mapped sequences.

2. Graph constructed from long read or sequence data, extension with long reads or sequences
--------------------------------------------------------------------------------------------

Here, our recommendation is to actually rebuild the graph with `PGGB <https://github.com/pangenome/pggb>`_. One could
use `Graphaligner <https://github.com/maickrau/GraphAligner>`_ to align the long
sequences against the graph and then use ``vg augment`` to extend the already existing graph, but that would be
comparatively inexact and the resolutions of complex regions might drop dramatically.
A reference-biased method would be `Minigraph <https://github.com/lh3/minigraph>`_ followed by `Cactus <https://github.com/glennhickey/progressiveCactus>`_.

3. Graph constructed from short read or sequence data, extension with short reads or sequences
----------------------------------------------------------------------------------------------

Map the sequences to the graph with ``vg giraffe``. The resulting `GAM <https://github.com/vgteam/vg/wiki/File-Formats#gam-graph-alignment--map-vgs-bam>`_
file can be used with ``vg augment`` to extend the existing graph with the mapped sequences.

4. Graph constructed from short read or sequence data, extension with long reads or sequences
---------------------------------------------------------------------------------------------

Use `Graphaligner <https://github.com/maickrau/GraphAligner>`_ to align the long sequences against the graph and then use
``vg augment`` to extend the already existing graph.

All of the above methods produce a pangenome graph in GFAv1 format which can then be analysed with ``odgi``.

Why is ``odgi`` strictly limited to GFAv1? Why does it not support GFAv2 or rGFA?
=================================================================================

Although `GFAv2 is a superset of GFAv1 <http://gfa-spec.github.io/GFA-spec/GFA2.html#backward-compatibility-with-gfa-1>`_,
GFAv2 was specifically designed for assembly graphs. The fields required to losslessly represent a variation graph are
already specified in the more frequently used GFAv1.
The `rGFA <https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-reference-gfa-rgfa-format>`_ format requires a genomic sequence to be the reference
sequence upon which all other sequences are related to. In GFAv1 we don't have that limitation and this is fundamental
to implement reference-free approaches.

How is heterozygosity handled by ``odgi``? How polyploidy?
==========================================================

The GFA format doesn’t store the metadata information. To overcome this limit, we store biosample information in the
sequence names that become the path names in the graph, by following the `PanSN-spec convention <https://github.com/pangenome/PanSN-spec>`_.
In more detail, we apply the following sequence naming scheme for sequences:

``[sample_name][delimiter][haplotype_id][delimiter][contig_or_scaffold_name]``

Where each field is optional. For instance, by using the character ‘#’ as delimiter, the sequence name 'HG002#1#ctg1234’
names ‘ctg1234’ on the first haplotype (or phase group) of the HG002 individual, while ‘HG002#2#ctg9876’ is contig ‘ctg9876’
on the other haplotype of the same individual. This can be naturally extended and applied for polyploid species as well.
To give a concrete example: If one only wants to work with a graph containing the associated haplotypes,
:ref:`odgi extract` can be restricted to the desired haplotypes with the `-p[FILE],--paths-to-extract=[FILE]` parameter.

How does :ref:`odgi position`'s GFF liftover work?
==================================================

The GFF file contains annotations for one or more paths in the graph. For each annotation, we know the start and end
within that path. So we can annotate all nodes that are visited by such a path range with the information from the
`attribute` field. If there are overlapping features, we append the annotation for each node. Using the same coloring
schema as in :ref:`odgi viz` we generate a color for each annotated node by its collected annotation.

If a subgraph was as a result from e.g. :ref:`odgi extract`, the path names are usually in the form of
``name:start-end``. :ref:`odgi position` is able to automatically detect this and adjust the positions given in the GFF
on the fly to the new positions given in the subgraph. For each GFF entry, it just subtracts the “missing” number of
nucleotides from the `start` and end `field`. That’s how we adjust for the subgraph annotation.


