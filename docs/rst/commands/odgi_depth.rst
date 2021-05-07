.. _odgi depth:

#########
odgi depth
#########

find the depth of graph as defined by query criteria

SYNOPSIS
========

**odgi depth** [**-i, –input**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi depth(1) command finds the depth of graph as defined by query
criteria.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –input**\ =\ *FILE*
| Compute path depths in this graph.

Depth Options
-------------

| **-s, –subset-paths**\ =\ *FILE*
| Compute the depth considering only the paths specified in the FILE;
  the file must contain one path name per line and a subset of all paths
  can be specified.

| **-r, –path**\ =\ *STRING*
| Compute the depth of the given path in the graph.

| **-R, –paths**\ =\ *FILE*
| Compute depth for the paths listed in FILE.

| **-g, –graph-pos**\ =\ *[[node_id][,offset[,(+|-)]\ *\ **]**\ *]*
| Compute the depth at the given node, e.g. 7 or 3,4 or 42,10,+ or
  302,0,-.

| **-G, –graph-pos-file**\ =\ *FILE*
| A file with graph path position per line.

| **-p, –path-pos**\ =\ *[[path_name][,offset[,(+|-)]\ *\ **]**\ *]*
| Return depth at the given path position e.g. chrQ or chr3,42 or
  chr8,1337,+ or chrZ,3929,-.

| **-F, –path-pos-file**\ =\ *FILE*
| A file with one path position per line.

| **-b, –bed-input**\ =\ *FILE*
| A BED file of ranges in paths in the graph.

| **-d, –graph-depth**
| Compute the depth on each node in the graph.

| **-d, –search-radius**\ =\ *STRING*
| Limit coordinate conversion breadth-first search up to DISTANCE bp
  from each given position [default: 10000].

Threading
---------

| **-t, –threads**\ =\ *N*
| Number of threads to use.

Program Information
-------------------

| **-h, –help**
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
