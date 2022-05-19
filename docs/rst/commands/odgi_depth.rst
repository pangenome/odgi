.. _odgi depth:

#########
odgi depth
#########

Find the depth of a graph as defined by query criteria. Depth can be node-depth computed by counting all defined paths that run through a node. By default it prints in a tab-delimited format path, start, end, and mean.depth to stdout.

SYNOPSIS
========

**odgi depth** [**-i, --input**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi depth command finds the depth of graph as defined by query
criteria. Without specifying any non-mandatory options, it prints in a tab-delimited
format *path*, *start*, *end*, and *mean.depth* to stdout.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --input**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

Depth Options
-------------

| **-s, --subset-paths**\ =\ *FILE*
| Compute the depth considering only paths specified in *FILE*.
  The file must contain one path name per line and a subset of all paths
  can be specified. If a step is of a path of the given list, it is taken into account when calculating a node's depth.

| **-r, --path**\ =\ *STRING*
| Compute the depth of the given path *STRING* in the graph.

| **-R, --paths**\ =\ *FILE*
| Report the depth only for the paths listed in *FILE*.

| **-g, --graph-pos**\ =\ *[[node_id][,offset[,(+|-)]\ *\ **]**\ *]*
| Compute the depth at the given node e.g. 7 or 3,4 or 42,10,+ or
  302,0,-.

| **-G, --graph-pos-file**\ =\ *FILE*
| A file with one graph path position per line.

| **-p, --path-pos**\ =\ *[[path_name][,offset[,(+|-)]\ *\ **]**\ *]*
| Return depth at the given path position e.g. chrQ or chr3,42 or
  chr8,1337,+ or chrZ,3929,-.

| **-F, --path-pos-file**\ =\ *FILE*
| A file with one path position per line.

| **-b, --bed-input**\ =\ *FILE*
| A BED file of ranges in paths in the graph.

| **-d, --graph-depth-table**
| Compute the depth and unique depth on each node in the graph, writing a table by node:
*node.id*, *depth*, and *depth.uniq*.

| **-v, --graph-depth-vec**
| Compute the depth on each node in the graph, writing a vector by base in one line.

| **-D, --path-depth**
| Compute a vector of depth on each base of each path. Each line consists of a path name
 and subsequently the space-separated depth of each base.

| **-a, --self-depth**
| Compute the depth of the path versus itself on each base in each path. Each line consists of a path name
 and subsequently the space-separated depth of each base.

| **-S, --summarize**
| Provide a summary of the depth distribution in the graph. In a tab-delimited format it
 prints to stdout: *node.count*, *graph.length*, *step.count*, *path.length*,
  *mean.node.depth* (step.count/node.count), and *mean.graph.depth* (path.length/graph.length).

| **-w, --windows-in**\ =\ *LEN:MIN:MAX*
| Print to stdout a BED file of path intervals where the depth is between *MIN* and
 *MAX*, merging the ranges not separated by more then *LEN* bp.

| **-W, --windows-out**\ =\ *LEN:MIN:MAX*
| Print to stdout a BED file of path intervals where the depth is outside *MIN* and
 *MAX*, merging the ranges not separated by more then *LEN* bp.

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use in parallel operations.

Processing Information
----------------------

| **-P, --progress**
| Print information about the operations and the progress to stderr.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi depth**.

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
