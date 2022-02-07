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


