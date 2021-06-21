.. _tools_using_odgi:

#################
Tools using ODGI
#################

The Pangenome Graph Builder (PGGB)
==================================

`PGGB <https://github.com/pangenome/pggb>`_ renders a collection of sequences into a pangenome graph. Its goal is to build a graph that is locally directed and
acyclic while preserving large-scale variation. Maintaining local linearity is important for the interpretation,
visualization, and reuse of pangenome variation graphs.
It uses three phases:

  1. `wfmash <https://github.com/ekg/wfmash>`_: The probabilistic mash-map mapper ``wfmash`` is used to create all versus all alignments of all input segments.

  2. `seqwish <https://github.com/ekg/seqwish>`_: The pangenome graph is induced from the alignments and emitted as a GFAv1.

  3. `smoothxg <https://github.com/pangenome/smoothxg>`_: The graph is then sorted with a form of multi-dimensional scaling in 1D, groomed, and topologically ordered locally. The 1D order is then broken into "blocks" which are "smoothed" using the partial order alignment (POA) algorithm implemented in abPOA.

The linearizaton and simplification of a variation graph is possible using the various sorting algorithms in :ref:`odgi sort`.
The final pangenome graph is translated to the ODGI format using :ref:`odgi build`. :ref:`odgi stats` reports statistics
of the created graph and a 1D visualization from :ref:`odgi viz` allows diagnostic analyses. The 2D layouts and drawings by
:ref:`odgi layout` and :ref:`odgi draw` complement the visual analytic tools.

Burning complex regions of a pangenome
======================================

Generating a pangenome graph from over 40 human assemblies is a challenging task. When using ``PGGB``, centromers appear
as complex regions in the resulting graph. For the ease of downstream analysis, such regions can be burned away. An example
of such a script can be found at `burn_extract.sh <https://github.com/pangenome/HPRCy1v2/blob/main/burn_extract.sh>`_ .

:ref:`odgi depth` identifies low and high coverage regions which can be removed from the graph with :ref:`odgi extract`.
:ref:`odgi sort` ensures local and global linearity. :ref:`odgi stats`, :ref:`odgi depth`, and :ref:`odgi degree` deliver
in detail statistics about the burned graph. A visual analysis is possible with :ref:`odgi viz`. To go full circle,
:ref:`odgi view` projects the pangenome graph back to GFAv1 format.

Nextflow nf-core/pangenome pipeline
===================================

A Nextflow version of ``PGGB`` is currently developed on nf-core. The `nf-core/pangenome <https://github.com/nf-core/pangenome>`_ pipeline makes use of the same
``odgi`` tools as the original ``PGGB``.

Interactive 2D visualization of very large graphs - gfaestus
============================================================

`gfaestus <https://github.com/chfi/gfaestus>`_ is a very potent 2D interactive visualization tool for variation graphs.
It depends on the 2D TSV layout produced by :ref:`odgi layout`. For further information take a look at the :ref:`sorting-layouting` section.
