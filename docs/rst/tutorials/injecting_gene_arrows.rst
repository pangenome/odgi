.. _untangling_the_pangenome:

######################
Injecting gene arrows
######################

========
Synopsis
========

A pangenome graph represents the alignment of many genome sequences.
By embedding gene annotations into the graph as `paths`, we "align" them with all other paths.

We call this embedding `injection`.
It's implemented in ``odgi inject``.
Injection makes new paths that are equal to ranges (or alignments) on graph paths.
Of course, these trivially match their source path.
For instance, if we take a gene annotation (e.g. in BED format) on GRCh38, then "inject" it into the graph, we'll get a new path that exactly matches a corresponding sub-path of GRCh38.

However, things get interesting when we use these new gene paths as reference targets for ``odgi untangle``.
Untangling shows us how each `query` path maps onto a set of `targets`.
Using our genes as targets will show the order and relative orientation of each path in the graph relative to a set of gene annotations.
This process projects our gene annotations onto overlapping genomes.

In this tutorial, we demonstrate injection and untangling on the human ``C4`` locus.
Our goal is to make gene arrow maps to show the different haplotypes relative to copies of the C4 gene in the GRCh38 reference (C4A and C4B).
We use the `gggenes R package <https://github.com/wilkox/gggenes>`_ to plot gene orders across haplotypes in this locus.

=====
Steps
=====

-----------------------------
Build the C4 graph
-----------------------------

Assuming that your current working directory is the root of the ``odgi`` project, to construct an ``odgi`` file from the
``C4`` dataset in ``GFA`` format, execute:

.. code-block:: bash

    odgi build -g test/chr6.C4.gfa -o chr6.C4.og

The command creates a file called ``chr6.C4.og``, which contains the input graph in ``odgi`` format. This graph contains
90 contigs from 88 haploid human genome assemblies from 44 individuals plus the grch38 and chm13 reference sequences.
The contigs cover the `Complement Component 4 (C4) region <http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6%3A31972057%2D32055418&hgsid=1211002763_taymHTHRXlpfFiqu51J8nWGyKv67>`_.
In humans, C4 is a protein involved in the intricate complement system, originating from the human leukocyte antigen (HLA) system.


-----------------------------------
Injecting genes
-----------------------------------

We start with gene annotations against the GRCh38 reference.
Our annotations are against the full `grch38#chr6`, in ``test/chr6.C4.bed``:

.. code-block:: none

   grch38#chr6	31982057	32002681	C4A
   grch38#chr6	32014795	32035418	C4B
   grch38#chr6	31985101	31991058	C4A_HERV-K
   grch38#chr6	32017839	32023796	C4B_HERV-K

However, the C4 locus graph ``chr6.c$.gfa`` is over the reference range ``grch38#chr6:31972046-32055647``.
``odgi extract`` uses `PanSN sequencing naming format <https://github.com/pangenome/PanSN-spec>`_, and so we can see the

.. code-block:: bash

   odgi paths -i chr6.C4.gfa -L | grep grc
   # grch38#chr6:31972046-32055647

We must adjust the annotations to match the subgraph to ensure that path names and coordinates exactly correspond between the BED and GFA.
We subtract ``31972046`` from both coordinates and adjust the path name to match that in the subgraph, yielding:

.. code-block:: bash

   <chr6.C4.bed awk '{ i=31972046; j=32055647; print $1":"i"-"j, $2-i, $3-i, $4 }' | tr ' ' '\t' >chr6.C4.adj.bed

Now, we can inject these into the graph:

.. code-block:: bash

   odgi inject -i chr6.C4.gfa -b chr6.C4.adj.bed -o - | odgi paths -i - -L | tail -4 | sort -V >chr6.C4.gene.names.txt

This shows that the annotations have been added as paths (``cat chr6.C4.gene.names.txt``):

.. code-block:: none

   C4A
   C4A_HERV-K
   C4B
   C4B_HERV-K

We can always pipe the output of ``odgi`` subcommands to each other, but in this case it will simplify things to save the graph with the injected gene paths:

.. code-block:: bash

   odgi inject -i chr6.C4.gfa -b chr6.C4.adj.bed -o chr6.C4.genes.og


-----------------------------
Visualize the C4 graph
-----------------------------

To visualize a subset of the graph, execute:

.. code-block:: bash

    # Select haplotypes
    odgi paths -i chr6.C4.og -L | grep 'chr6\|HG00438\|HG0107\|HG01952\|C4' > chr6.C4.selected_paths.txt

    odgi viz -i chr6.C4.og -o chr6.C4.selected_paths.png -c 12 -w 100 -y 50 -p chr6.C4.selected_paths.txt -m -B Spectral:4

To obtain the following PNG image:

.. image:: /img/chr6.C4.selected_paths_inject.png

The selected paths (2 reference genomes and 6 haplotypes of 3 individuals) are colored by path depth.
We additionally see the C4 annotation paths at the bottom of the visualization.
Several color palettes are available (see :ref:`odgi viz` documentation for more information), with the default ``Spectral`` palette suitable for examining collapsed repeats in the graph.
(Here we use the ``Spectral:4``-color version to increase readability, but ``Spectral:11`` is default with ``odgi viz -m``.)

Human C4 exists as 2 functionally distinct genes, C4A and C4B, which both vary in structure and copy number (`Sekar et al., 2016 <https://doi.org/10.1038/nature16549>`_).
By injecting annotations into the graph, we can see where these copies fit (bottom 4 path rows in the image above).
The longer link on the bottom indicates that the copy number status varies across the haplotypes represented in the pangenome.
Moreover, C4A and C4B genes segregate in both long and short genomic forms, distinguished by the presence or absence of a
human endogenous retroviral (HERV) sequence, as also highlighted by the short nested link on the left of the image.

Coloring by path depth, we can see that the two references present two different allele copies of the C4 genes,
both of them including the HERV sequence. The entirely grey paths have one copy of these genes. ``HG01071#2`` presents 3 copies of the locus (orange),
of which one contains the HERV sequence (gray in the middle of the orange). In ``HG01952#1``, the HERV sequence is absent.

-----------------------------------
Untangling to obtain a gene arrow map
-----------------------------------

We now use the gene names and the ``gggenes`` output format from ``odgi untangle`` to obtain a gene arrow map!
We use `-j 0.5` to filter out low-quality matches.

.. code-block:: bash

   odgi untangle -R chr6.C4.gene.names.txt -i chr6.C4.genes.og -j 0.5 -t 16 | grep '^mol\|HG00438#2\|HG0107\|HG01952#1\|chm13' >chr6.C4.gene.gggenes.tsv


Take a look at the rows in the ``chr6.C4.untangle.bed`` file for the ``HG00438`` and ```HG01071`` individuals:

.. code-block:: bash

    cat <(head chr6.C4.untangle.bed -n 1) <(grep 'HG00438\|HG01071' chr6.C4.untangle.bed) | column -t

.. code-block:: none




For each segment in the query (``query.name``, ``query.start``, and ``query.end`` columns), the best match on the reference is reported
(``ref.name``, ``ref.start``, and ``ref.end``), with information about the quality of the match (``score``), the strand (``inv``),
the copy number status (``self.cov``), and its rank over all possible matches (``n.th``).

To obtain a visualization of the output for the ``HG00438`` and ```HG01071`` individuals, execute:

.. code-block:: R

    library(ggplot2)
    x <- read.table('chr6.C4.untangle.bed', sep = '\t', header = T, comment.char="$")
    x$query.name <- gsub(":.*","",x$query.name)
    x$query.name <- gsub("#J.*","",x$query.name)

    ggplot(
      subset(x, query.name %in% c(
        "grch38#chr6",
        "HG00438#1",
        "HG00438#2",
        "HG01071#1",
        "HG01071#2"
        )
      ), aes(x=query.start, xend=query.end, y=ref.start, yend=ref.end)) +
        geom_segment(size=0.3) +
        facet_grid(. ~ query.name) +
        coord_fixed() +
        theme(
          text = element_text(size = 12.6),
          axis.text.x = element_text(size = 12, angle = 90),
          axis.text.y = element_text(size = 12),
        )  +
          xlab("Query start") +
          ylab("Reference start")

    ggsave('chr6.C4.untangle.png', width = 32, height = 8,  units = "cm", dpi = 300,  bg = "transparent")


To obtain the following PNG image:

.. image:: /img/chr6.C4.untangle.png

The plots show the copy number status of the haplotypes in the C4 region with respect to the grch38 reference sequence.
On the grch38 reference, `C4A precedes C4B, and both are in single copy <http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6%3A31972057%2D32055418&hgsid=1211002763_taymHTHRXlpfFiqu51J8nWGyKv67>`_.
``odgi untangle``'s output makes then clear, for example, that in ``HG00438`` the C4A gene is missing in both haplotypes, while ``HG01071#2``
has two copies of C4B.
