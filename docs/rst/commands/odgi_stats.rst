.. _odgi stats:

#########
odgi stats
#########

metrics describing variation graphs

SYNOPSIS
========

**odgi stats** [**-i, –idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi stats(1) command produces statistics of a variation graph.
Among other metrics, it can calculate the #nodes, #edges, #paths and the
total nucleotide length of the graph. Various histogram summary options
complement the tool. If [**-B, –bed-multicov**\ =\ *BED*] is set, the
metrics will be produced for the intervals specified in the BED.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to create statistics
  from. The file name usually ends with *.og*.

Summary Options
---------------

| **-S, –summarize**
| Summarize the graph properties and dimensions. Print to stdout the
  #nucleotides, #nodes, #edges and #paths of the graph.

| **-W, –weak-connected-components**
| Shows the properties of the weakly connected components.

| **-L, –self-loops**
| Number of nodes with a self-loop.

| **-N, –nondeterministic-edges**
| Show nondeterministic edges (those that extend to the same next base).

| **-b, –base-content**
| Describe the base content of the graph. Print to stdout the #A, #C, #G
  and #T of the graph.

Sorting goodness evaluation
---------------------------

| **-l, –mean-links-length**
| Calculate the mean links length. This metric is path-guided and
  computable in 1D and 2D.

| **-g, –no-gap-links**
| Don’t penalize gap links in the mean links length. A gap link is a
  link which connects two nodes that are consecutive in the linear
  pangenomic order. This option is specifiable only to compute the mean
  links length in 1D.

| **-s, –sum-path-nodes-distances**
| Calculate the sum of path nodes distances. This metric is path-guided
  and computable in 1D and 2D. For each path, it iterates from node to
  node, summing their distances, and normalizing by the path length. In
  1D, if a link goes back in the linearized viewpoint of the graph, this
  is penalized (adding 3 times its length in the sum).

| **-d, –penalize-different-orientation**
| If a link connects two nodes which have different orientations, this
  is penalized (adding 2 times its length in the sum).

| **-c, –coords-in**
| File containing the layout coordinates of the succinct variation graph
  specified as input. The file name usually ends with *.lay*. When the
  layout coordinates are provided, the mean links length and the sum
  path nodes distances statistics are evaluated in 2D, else in 1D.

| **-P, –path-statistics**
| Display the statistics (mean links length or sum path nodes distances)
  for each path.

Threading
---------

| **-t, –threads**\ =\ *N*
| Number of threads to use.

Program Information
-------------------

| **-h, –help**
| Print a help message for **odgi stats**.

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
