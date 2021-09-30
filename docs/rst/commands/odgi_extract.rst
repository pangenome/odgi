.. _odgi extract:

#########
odgi extract
#########

Extract subgraphs or parts of a graph defined by query criteria.

SYNOPSIS
========

**odgi extract** [**-f, --input-graphs**\ =\ *FILE*] [**-o,
--out**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi extract command extracts parts of the graph defined by
query criteria.

OPTIONS
=======

MANDATORY OPTIONS
-----------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

Graph Files IO
--------------

| **-o, --out**\ =\ *FILE*
| Store all subgraphs in this FILE. The file name usually ends with
  *.og*.

Extract Options
---------------

| **-s, --split-subgraphs**
| Instead of writing the target subgraphs into a single graph, write one
  subgraph per given target to a separate file named
  ``path:start-end.og`` (0-based coordinates).

| **-I, --inverse**
| Extract the parts of the graph that do not meet the query criteria.

| **-l, --node-list**\ =\ *FILE*
| A file with one node id per line. The node specified will be extracted
  from the input graph.

| **-n, --node**\ =\ *ID*
| Specify a single node ID from which to begin our traversal.

| **-c, --context-steps**\ =\ *N*
| The number of steps (nodes) *N* away from our initial subgraph that we should
  collect.

| **-L, --context-bases**
| The number of bases *N* away from our initial subgraph that we should collect.

| **-r, --path-range**
| Find the node(s) in the specified path range TARGET=path[:pos1[-pos2]]
  (0-based coordinates).

| **-r, --bed-file**\ =\ *FILE*
| Find the node(s) in the path range(s) specified in the given BED FILE

| **-E, --full-range**
| Collects all nodes in the sorted order of the graph in the min and max
  position touched by the given path ranges. This ensures that all the paths of the subgraph are not split by node, but that the nodes are laced together again. Comparable to **-R, --lace-paths=FILE**, but specifically for all paths in the resulting subgraph. Be careful to use it with
  very complex graphs.

| **-p, --paths-to-extract**\ =\ *FILE*
| List of paths to consider in the extraction. The *FILE* must contain one
  path name per line and a subset of all paths can be specified.

| **-R, --lace-paths**\ =\ *FILE*
| List of paths to fully retain in the extracted graph. Must contain one
  path name per line and a subset of all paths can be specified.

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
| Print a help message for **odgi extract**.

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

