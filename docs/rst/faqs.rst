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

Why are even large pangenome graphs expected to be sparse?
==========================================================

“A dense graph is a graph in which the number of edges is close to the maximum number of edges.” (https://en.wikipedia.org/wiki/Dense_graph).
Consequently, a sparse graph is a graph in which the number of edges is much less than the possible number of edges.
As we allow self-edges and we have a bidirected graph, the number of maximum edges can be calculated with

.. math::

   2* {\sum_{i=1}^{n}i-1} + n

where *n* is the number of nodes in the graph. One would classify a graph as sparse if the number of edges is at most half
the number of maximum edges of that graph. The HTTexon1 graph from Figure 3 of the paper has 35 nodes and 56 edges.
So the maximum number of edges is 1225. The graph clearly is sparse. The centromere of a 90 haplotype human pangenome
chromosome 8 (https://www.nature.com/articles/s41586-021-03420-7) graph has 377123 nodes and 560986 edges. The maximum
number of edges is 142,221,757,129. The graph clearly is very sparse. Collaborators are currently building a Cannabis sativa
pangenome graph from 12 haplotypes. Their current chromosome 7 graph has ~2M nodes and ~2.8M edges. The maximum number
of edges would be 40,000,000,000. The graph clearly is very sparse. In the evaluation of the complex graphs above, we
only observe very sparse graphs. Even the C. sativa graph, although each genome is to be expected to consist of ~75%
of repetitive elements (https://www.nature.com/articles/s41438-020-0295-3), is very sparse.

Genomic obesity (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC157029/) leads to more repetitive elements in a graph. In a variation graph, a repetitive sequence is usually
compressed into one node with some edges connecting to the other nodes. The graph would get dense if a repetitive node
would be connected to several million other nodes in the graph. We currently can’t think of any biological data that
would lead to such a phenomenon. Therefore,  we argue that a large number of repetitive elements does not lead to a
significantly higher number of edges, but to a significantly higher number of steps in the graph. And that’s one of
the major use cases ``odgi`` was designed for: its Step structure and parallel path processing capabilities allow it to
work with a very large number of paths and steps in a graph.

At some point in time, adding more and more genomes to an already very large pangenome will lead to (i) a core pangenome
that will basically never change, and (ii) some individual genomic sequences that will add more sequence to the pangenome.
But, also taking the arguments of the paragraph above into account, adding new sequences even from complex regions like
centromeres into a large pangenome graph won’t lead to a dense, but a sparse graph, too. Ultimately, the construction
method and the variation encoded in a pangenome graph have the greatest influence on the sparseness of a graph.
Clearly, we can make a graph consisting of a very small number of nodes that represent, e.g., all extant 5-mers.
This graph will not be sparse, but it will also be very different and serve a different research objective than
pangenome, assembly, and multiple alignment graphs typically used in the research community.

Can ODGI accurately represent repeats?
======================================

Yes! The transformation of a graph in GFAv1 format to ``odgi``’s binary format is lossless. Indeed, a graph
in ``odgi`` format fully represents all nodes, edges, paths, sequences, and any kind of variation present in the input GFAv1
file. Of note, the initial graph construction method itself, for example PGGB (https://github.com/pangenome/pggb) or
minigraph (Li et al., 2020), determines the encoding of the repeats in the input graph.

How does the sorting and change of the node order work in general?
==================================================================

The sort order of the graph is the order in which nodes are enumerated. We can assign new node IDs to change the sort
order. We find that, because they are typically very sparse, sorting pangenome graphs can help to reveal underlying
structures and patterns of variation. This is key for visualization and interpretation.

Most subcommands in ``odgi`` require and verify that the input graph’s node identifiers (IDs) are optimized, that is
compacted from 1 to *N* where *N* is the number of nodes in the graph. If this assumption is violated, :ref:`odgi sort` provides
functionality to optimize the graph. This means that the first node identifier (ID) starts at 1 and the last node ID is
the number of nodes. All sorting operations update the graph in place with an efficient ID rewriting algorithm.
The graph is then updated in place. First, the node identifiers are normalized (from 1 to number of nodes) including
the adjustment of the edges. Second, path information, including both path metadata that points into the start and end
steps of the path, plus each step of every path, is updated, too. We point out that changing the node order does not
change our coordinate systems based on paths. These will now refer to a new node ordering.

Does ODGI groom remove true minor variants?
===========================================

No, minor variants would not be removed. The grooming process is lossless with respect to graph content and overall topology:
what is altered is the local orientation of the assemblies in the pangenome graph, with the aim of simplifying the graph
structure for easier downstream analyses. Grooming works to simplify the representation of inversions, to require fewer
edges that go between the two strands of the graph.
