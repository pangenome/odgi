.. _presence_absence_variants:

####################################
Identify presence/absence variants
####################################

========
Synopsis
========

The term presence/absence variation (PAV) is used to describe sequences that are present in one genome, but
entirely missing in another genome, and is an important source of genetic divergence and diversity.

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

-----------------------
Presence/absence variants (PAVs)
-----------------------

Any path in the graph can be used as a reference to identify PAVs. In this example, we have chosen the ``chm13__LPA__tig00000001``
path. To obtain contains 1000 bp interval windows across the chosen reference, execute:

.. code-block:: bash

    odgi paths -i LPA.og -f | grep 'chm13__LPA__tig00000001' -A 1 > reference.fa
    samtools faidx reference.fa
    bedtools makewindows -g <(cut -f 1,2 reference.fa.fai) -w 1000 > LPA.w1kbp.bed

To identify the PAVs, execute:

.. code-block:: bash

    odgi pav -i LPA.og -b LPA.w1kbp.bed > LPA.w1kbp.pavs.txt

By default, ``odgi pav`` prints to stdout a matrix with the PAVs ratios.
To take a look at the first rows and columns of the matrix in the ``LPA.w1kbp.pavs.txt`` file, execute:

.. code-block:: bash

    head LPA.w1kbp.pavs.txt | cut -f 1-8 | column -t

.. code-block:: none

    chrom                    start  end   name  chm13__LPA__tig00000001  HG002__LPA__tig00000001  HG002__LPA__tig00000005  HG00733__LPA__tig00000001
    chm13__LPA__tig00000001  0      1000  .     1                        0                        0                        0
    chm13__LPA__tig00000001  1000   2000  .     1                        0                        0                        0
    chm13__LPA__tig00000001  2000   3000  .     1                        0                        0                        0
    chm13__LPA__tig00000001  3000   4000  .     1                        0                        0                        0
    chm13__LPA__tig00000001  4000   5000  .     1                        0                        0                        0
    chm13__LPA__tig00000001  5000   6000  .     1                        0.4156                   0.91101                  0.00091743
    chm13__LPA__tig00000001  6000   7000  .     1                        1                        1                        0.80339
    chm13__LPA__tig00000001  7000   8000  .     1                        0.99811                  0.99906                  0.98491
    chm13__LPA__tig00000001  8000   9000  .     1                        1                        1                        0.99466


The ``chrom``, ``start``, ``end``, and ``name`` columns are filled with the values in the corresponding columns in the
input ``BED`` format file. In this example, the region ``chm13__LPA__tig00000001:5000-6000`` is covered at ``41.56%`` in the
``HG002__LPA__tig00000001`` contig and it is almost absent in ``HG00733__LPA__tig00000001``.

To emit a binary PAV matrix, execute:

.. code-block:: bash

    odgi pav -i LPA.og -b LPA.w1kbp.bed -B 0.5 > LPA.w1kbp.pavs.binary.txt
    head LPA.w1kbp.pavs.binary.txt | cut -f 1-8 | column -t

.. code-block:: none

    chrom                    start  end   name  chm13__LPA__tig00000001  HG002__LPA__tig00000001  HG002__LPA__tig00000005  HG00733__LPA__tig00000001
    chm13__LPA__tig00000001  0      1000  .     1                        0                        0                        0
    chm13__LPA__tig00000001  1000   2000  .     1                        0                        0                        0
    chm13__LPA__tig00000001  2000   3000  .     1                        0                        0                        0
    chm13__LPA__tig00000001  3000   4000  .     1                        0                        0                        0
    chm13__LPA__tig00000001  4000   5000  .     1                        0                        0                        0
    chm13__LPA__tig00000001  5000   6000  .     1                        0                        1                        0
    chm13__LPA__tig00000001  6000   7000  .     1                        1                        1                        1
    chm13__LPA__tig00000001  7000   8000  .     1                        1                        1                        1
    chm13__LPA__tig00000001  8000   9000  .     1                        1                        1                        1

With ``B`` is specified to emit a binary matrix, with 1 if the PAV ratio is greater than or equal to the specified
threshold (``0.5`` in the example), else 0.

If needed, it is possible to group paths. For this, we need to prepare a file that specifies for each path the group it
belongs to. In the ``LPA`` pangenome graph, the first part of each path name indicates the sample name. Therefore, to
prepare such a file, execute:

.. code-block:: bash

    odgi paths -i LPA.og -L > LPA.paths.txt
    cut -f 1 -d '_' LPA.paths.txt > LPA.samples.txt
    paste LPA.paths.txt LPA.samples.txt > LPA.path_and_sample.txt

    head LPA.path_and_sample.txt -n 5 | column -t

.. code-block:: none

    chm13__LPA__tig00000001    chm13
    HG002__LPA__tig00000001    HG002
    HG002__LPA__tig00000005    HG002
    HG00733__LPA__tig00000001  HG00733
    HG00733__LPA__tig00000008  HG00733

Then, to group the PAVs by sample, execute:

.. code-block:: bash

    odgi pav -i LPA.og -b LPA.w1kbp.bed -B 0.5 -p LPA.path_and_sample.txt > LPA.w1kbp.pavs.binary.grouped_by_sample.txt

    head LPA.w1kbp.pavs.binary.grouped_by_sample.txt | column -t

.. code-block:: none

    chrom                    start  end   name  HG002  HG00733  HG01358  HG02572  NA19239  NA19240  chm13
    chm13__LPA__tig00000001  0      1000  .     0      0        0        1        0        0        1
    chm13__LPA__tig00000001  1000   2000  .     0      0        0        1        0        0        1
    chm13__LPA__tig00000001  2000   3000  .     0      0        0        1        0        0        1
    chm13__LPA__tig00000001  3000   4000  .     0      0        0        1        0        0        1
    chm13__LPA__tig00000001  4000   5000  .     0      0        0        1        0        0        1
    chm13__LPA__tig00000001  5000   6000  .     1      0        0        1        0        0        1
    chm13__LPA__tig00000001  6000   7000  .     1      1        1        1        0        0        1
    chm13__LPA__tig00000001  7000   8000  .     1      1        1        1        1        0        1
    chm13__LPA__tig00000001  8000   9000  .     1      1        1        1        1        0        1
