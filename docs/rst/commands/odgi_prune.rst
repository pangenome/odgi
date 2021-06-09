.. _odgi prune:

#########
odgi prune
#########

Remove complex parts of the graph.

SYNOPSIS
========

**odgi prune** [**-i , --idx**\ =\ *FILE*] [**-o, --out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi prune command can remove complex parts of a graph. One can
drop paths, nodes by a certain kind of edge coverage, edges and graph
tips. Specifying a kmer length and a maximum number of furcations, the
graph can be broken at edges not fitting into these conditions.

OPTIONS
=======

MANDATORY OPTIONS
-------------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*.

| **-o, --out**\ =\ *FILE*
| Write the pruned graph in ODGI format to *FILE*. A file ending with *.og* is recommended.

Kmer Options
------------

| **-k, --kmer-length**\ =\ *N*
| The length of the kmers to consider.

| **-e, --max-furcations**\ =\ *N*
| Break at edges that would induce *N* many furcations in a kmer.

Node Options
------------

| **-d, --max-degree**\ =\ *N*
| Remove nodes that have a higher node degree than *N*.

| **-c, --min-coverage**\ =\ *N*
| Remove nodes covered by fewer than *N* number of path steps.

| **-C, --max-coverage**\ =\ *N*
| Remove nodes covered by more than *N* number of path steps.

| **-T, --cut-tips**
| Remove nodes which are graph tips.

Edge Options
------------

| **-E, --edge-coverage**
| Remove edges outside of the minimum and maximum coverage rather than
  nodes. Only set this argument in combination with [**-c,
  –min-coverage**\ =\ *N*] and [**-C, --max-coverage**\ =\ *N*].

| **-b, --best-edges**\ =\ *N*
| Only keep the *N* most covered inbound and output edges of each node.

Step Options
------------

| **-s, --expand-steps**\ =\ *N*
| Also include nodes within this many steps of a component passing the prune thresholds.

| **-l, --expand-length**\ =\ *N*
| Also include nodes within this graph nucleotide distance of a component passing the prune thresholds.

Path Options
------------

| **-p, --expand-path-length**\ =\ *N*
| Also include nodes within this path length of a component passing the prune thresholds.

| **-D, --drop-paths**
| Remove the paths from the graph.

| **-m, --cut-tips-min-depth**\ =\ *N*
| Remove nodes which are graph tips and have less than *N* path depth.

| **-I, --remove-isolated**
| Remove isolated nodes covered by a single path.

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for parallel operations.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi prune**.

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
