.. _odgi bin:

#########
odgi bin
#########

Binning of pangenome sequence and path information in the graph.

SYNOPSIS
========

**odgi bin** [**-i, --idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

| The odgi bin command bins a given variation graph. The pangenome
  sequence, the one-time traversal of all nodes from smallest to largest
  node identifier, can be summed up into bins of a specified size. For
  each bin, the path metainformation is summarized. This enables a
  summarized view of gigabase scale graphs. Each step of a path is a bin
  and connected to its next bin via a link. A link has a start bin
  identifier and an end bin identifier.
| The concept of odgi bin is also applied in :ref:`odgi viz`. A demonstration of how the odgi
  bin JSON output can be used for an interactive visualization is
  realized in the `Pantograph <https://graph-genome.github.io/>`__
  project. Per default, odgi bin writes the bins to stdout in a
  tab-delimited format: **path.name**, **path.prefix**, **path.suffix**,
  **bin** (bin identifier), **mean.cov** (mean coverage of the path in
  this bin), **mean.inv** (mean inversion rate of this path in this
  bin), **mean.pos** (mean nucleotide position of this path in this
  bin), **first.nucl** (first nucleotide position of this path in this
  bin), **last.nucl** (last nucleotide position of this path in this
  bin). These nucleotide ranges might span positions that are not
  present in the bin. Example: A range of 1-100 means that the first
  nucleotide has position 1 and the last has position 100, but
  nucleotide 45 could be located in another bin. For an exact positional
  output, please specify [**-j, --json**].
| Running odgi bin in
  `HaploBlocker <https://github.com/tpook92/HaploBlocker>`__ mode, only
  arguments [**-b, --haplo-blocker**], [**-p[N],
  -–haplo-blocker-min-paths[N]**], and [**-c[N],
  -–haplo-blocker-min-coverage[N]**] are required. A TSV is printed to
  stdout: Each row corresponds to a node. Each column corresponds to a
  path. Each value is the coverage of a specific node of a specific
  path.

OPTIONS
=======

MANDATORY OPTIONS
-----------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

Bin Options
-----------

| **-n, --number-bins**\ =\ *N*
| The number of bins the pangenome sequence should be chopped up to.

| **-w, --bin-width**\ =\ *N*
| The bin width specifies the size of each bin.

| **-D, --path-delim**\ =\ *STRING*
| Annotate rows by prefix and suffix of this delimiter.

| **-a, --aggregate-delim**
| Aggregate on path prefix delimiter. Argument depends on [**-D,
  –path-delim**\ =\ *STRING*].

| **-j, --json**
| Print bins and links to stdout in pseudo JSON format. Each line is a
  valid JSON object, but the whole file is not a valid JSON! First, each
  bin including its pangenome sequence is printed to stdout per line.
  Second, for each path in the graph, its traversed bins including
  metainformation: **bin** (bin identifier), **mean.cov** (mean coverage
  of the path in this bin), **mean.inv** (mean inversion rate of this
  path in this bin), **mean.pos** (mean nucleotide position of this path
  in this bin), and an array of ranges determining the nucleotide
  position of the path in this bin. Switching first and last nucleotide
  in a range represents a complement reverse orientation of that
  particular sequence.

| **-s, --no-seqs**
| If [**-j, --json**] is set, no nucleotide sequences will be printed to
  stdout in order to save disk space.

| **-g, --no-gap-links**
| Don't include gap links in the output. We divide links into 2 classes:

1. The links which help to follow complex variations. They need to be
   drawn, else one could not follow the sequence of a path.

2. The links helping to follow simple variations. These links are called
   **gap-links**. Such links solely connecting a path from left to right
   may not be relevant to understand a path’s traversal through the
   bins. Therefore, when this option is set, the gap-links are left out
   saving disk space.

HaploBlocker Options
--------------------

| **-b, --haplo-blocker**
| Write a TSV to stdout formatted in a way ready for HaploBlocker: Each
  row corresponds to a node. Each column corresponds to a path. Each
  value is the coverage of a specific node of a specific path.

| **-p[N], --haplo-blocker-min-paths[N]**
| Specify the minimum number of paths that need to be present in the bin
  to actually report that bin. The default value is 1.

| **-c[N], --haplo-blocker-min-coverage[N]**
| Specify the minimum coverage a path needs to have in a bin to actually
  report that bin. The default value is 1.

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for parallel operations.

Processing Information
----------------------

| **-P, --progress**
| Print information about the operations and the progress to stderr.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi bin**.

..
	EXIT STATUS
	===========

	| **0**
	| Success.

	| **1**
	| Failure (syntax or usage error; parameter error; file processing
	  failure; unexpected error).

	BUGS
	====

	Refer to the **odgi** issue tracker at
	https://github.com/pangenome/odgi/issues.
