.. _odgi degree:

#########
odgi degree
#########

Describe the graph in terms of node degree.

SYNOPSIS
========

**odgi degree** [**-i, --idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi degree command describes the graph in terms of node degree.
In summarization mode, it shows the *node.count*, *edge.count*, *avg.degree*,
*min.degree*, and *max.degree*. One can also specify degree ranges streaming these into
a BED file.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

Summary Options
---------------

| **-s, --subset-paths**\ =\ *FILE*
| Compute the degree considering only the paths specified in the *FILE*.
  The file must contain one path name per line and a subset of all paths
  can be specified. If a step is of a path of the given list, it is taken into account when calculating a node's depth. Else not.

| **-r, --path**\ =\ *STRING*
| Compute the degree of the given path *STRING* in the graph.

| **-R, --paths**\ =\ *FILE*
| Report the degree only for the paths listed in *FILE*.

| **-g, --graph-pos**\ =\ *[[node_id][,offset[,(+|-)]\ *\ **]**\ *]*
| Compute the degree at the given node, e.g. 7 or 3,4 or 42,10,+ or
  302,0,-.

| **-G, --graph-pos-file**\ =\ *FILE*
| A file with one graph path position per line.

| **-p, --path-pos**\ =\ *[[path_name][,offset[,(+|-)]\ *\ **]**\ *]*
| Return degree at the given path position e.g. chrQ or chr3,42 or
  chr8,1337,+ or chrZ,3929,-.

| **-F, --path-pos-file**\ =\ *FILE*
| A file with one path position per line.

| **-b, --bed-input**\ =\ *FILE*
| A BED file of ranges in paths in the graph.

| **-d, --graph-degree-table**
| Compute the degree and unique degree on each node in the graph, writing a table by node:
*node.id*, *degree*, and *degree.uniq*.

| **-v, --graph-degree-vec**
| Compute the degree on each node in the graph, writing a vector by base in one line.

| **-D, --path-degree**
| Compute a vector of degree on each base of each path. Each line consists of a path name
 and subsequently the space-separated degree of each base.

| **-a, --self-degree**
| Compute the degree of the path versus itself on each base in each path. Each line consists of a path name
 and subsequently the space-separated degree of each base.

| **-S, --summarize-graph-degree**
| Summarize the graph properties and dimensions. Print to stdout the
  node.id and the node.degree.

| **-w, --windows-in**\ =\ *LEN:MIN:MAX*
| Print to stdout a BED file of path intervals where the degree is between *MIN* and *MAX*, merging the ranges not separated by more then *LEN* bp.

| **-W, --windows-out**\ =\ *LEN:MIN:MAX*
| Print to stdout a BED file of path intervals where the degree is outside *MIN* and *MAX*, merging the ranges not separated by more then *LEN* bp.

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
| Print a help message for **odgi degree**.

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
