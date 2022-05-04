.. _odgi heaps:

#########
odgi heaps
#########

Path pangenome coverage permutations.

SYNOPSIS
========

**odgi heaps** [**-i, --idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

Extract matrix of path pangenome coverage permutations for power law regression.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

Heaps Options
---------------

| **-p, --path-groups**\ =\ *FILE*
| Group paths as described in two-column *FILE*, with columns `path.name` and `group.name`.

| **-S, --group-by-sample**
| Following `PanSN <https://github.com/pangenome/PanSN-spec>`_ naming (`sample#hap#ctg`), group by sample (1st field).

| **-H, --group-by-haplotype**
| Following `PanSN <https://github.com/pangenome/PanSN-spec>`_ naming (`sample#hap#ctg`), group by haplotype (2nd field).

| **-b, --bed-targets**\ =\ *FILE*
| BED file over path space of the graph, describing a subset of the graph to consider.

| **-n, --n-permutations**\ =\ *N*
| Number of permutations to run.

| **-d, --min-node-depth**\ =\ *N*
| Exclude nodes with less than this path depth (default: 0).

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for parallel operations.

Processing Information
----------------------

| **-P, --progress**
| Print information about the components and the progress to stderr.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi heaps**.

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
