.. _odgi extract:

#########
odgi extract
#########

extract parts of the graph as defined by query criteria

SYNOPSIS
========

**odgi extract** [**-f, –input-graphs**\ =\ *FILE*] [**-o,
–out**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi extract(1) command extracts parts of the graph as defined by
query criteria.

OPTIONS
=======

Graph Files IO
--------------

| **-f, –input-graphs**\ =\ *FILE*
| File containing the succinct variation graph. The file name usually
  ends with *.og*.

| **-o, –out**\ =\ *FILE*
| Store all subgraph in this file. The file name usually ends with
  *.og*.

Extract Options
---------------

| **-s, –split-subgraphs**\ =\ *STRING*
| Instead of writing the target subgraphs into a single graph, write one
  subgraph per given target to a separate file named
  ``path:start-end.og`` (0-based coordinates).

| **-I, –inverse**
| Extract parts of the graph that do not meet the query criteria.

| **-l, –node-list**::_FILE\_
| A file with one node id per line. The node specified will be extracted
  from the input graph.

| **-n, –node**::_STRING\_
| A single node from which to begin our traversal.

| **-c, –context**::_NUMBER\_
| The number of steps away from our initial subgraph that we should
  collect.

| **-L, –use-length**
| Treat the context size as a length in bases (and not as a number of
  steps).

| **-r, –path-range**
| Find the node(s) in the specified path range TARGET=path[:pos1[-pos2]]
  (0-based coordinates)

| **-r, –bed-file**::_FILE\_
| Find the node(s) in the path range(s) specified in the given BED FILE

| **-E, –full-range**
| Collects all nodes in the sorted order of the graph in the min and max
  position touched by the given path ranges. Be careful to use it with
  very complex graphs.

| **-p, –paths-to-extract**::_FILE\_
| List of paths to consider in the extraction; the file must contain one
  path name per line and a subset of all paths can be specified.

| **-R, –lace-paths**::_FILE\_
| List of paths to fully retain in the extracted graph; must contain one
  path name per line and a subset of all paths can be specified.

Threading
---------

| **-t, –threads**\ =\ *N*
| Number of threads to use (to embed the subpaths in parallel).

Processing Information
----------------------

| **-P, –progress**
| Print information to stderr.

Program Information
-------------------

| **-h, –help**
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

