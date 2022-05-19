.. _odgi untangle:

#########
odgi untangle
#########

Project paths into reference-relative BEDPE, to decompose paralogy relationships.

SYNOPSIS
========

**odgi untangle** [**-i, --input**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi untangle command projects paths into a reference-relative BEDPE file, decomposing paralogy relationships. During this process, it is
capable of untangling loopy region resulting in linearized pairs of regions in the BEDPE file. A self dotplot assists in debugging and understanding
the untangle process.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --input**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

Untangling Options
----------------

| **-q, --query-path**\ =\ *NAME*
| Use this query path.

| **-r, --target-path**\ =\ *NAME*
| Use this target (reference) path.

| **-Q, --query-paths**\ =\ *FILE*
| Use query paths listed in *FILE*. Each line must contain one path.

| **-R, --target-paths**\ =\ *FILE*
| Use target (reference) paths list (one per line) in *FILE*.

| **-m, --merge-dist**\ =\ *N*
| Merge segments shorter than this length into previous segments.

| **-n, --n-best**\ =\ *N*
| Report up to the *N*th best target (reference) mapping for each query segment (default: *1*).

| **-j, --min-jaccard**\ =\ *N*
| Report target mappings >= the given jaccard threshold, with 0 <= *N* <= 1.0 (default: *0.0*).

| **-e, --cut-every**\ =\ *N*
| Cut every *N* base pairs of the sorted graph (default: *0/OFF*).

| **-p, --paf-output**
| Emit the output in PAF format.

| **-c, --cut-points-input**\ =\ *FILE*
A text file of node identifiers (one identifier per row) where to start the segment boundaries. When specified, no further starting points will be added.

| **-d, --cut-points-output**\ =\ *FILE*
Emit node identifiers where segment boundaries started (one identifier per row).

Debugging Options
-----------------

| **-s, --self-dotplot**
| Render a table showing the positional dotplot of the query against itself.

Step Index Options
------------------

| **-a, --step-index**\ =\ *FILE*
| Load the step index from this *FILE*. The file name usually ends with *.stpidx*. (default: build the step index from scratch with a sampling rate of 8).

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
| Print a help message for **odgi untangle**.

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
