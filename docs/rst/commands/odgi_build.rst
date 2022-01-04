.. _odgi build:

#########
odgi build
#########

Construct a dynamic succinct variation graph in ODGI format from a GFAv1

SYNOPSIS
========

**odgi build** [**-g, --gfa**\ =\ *FILE*] [**-o, --out**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi build command constructs a succinct variation graph from a
GFA. Currently, only GFAv1 is supported. For details of the format please
see https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md.

OPTIONS
=======

MANDATORY OPTIONS
-----------------

| **-g, --gfa**\ =\ *FILE*
| GFAv1 *FILE* containing the nodes, edges and paths to build a dynamic
  succinct variation graph from.

| **-o, --out**\ =\ *FILE*
| Write the dynamic succinct variation graph to this *FILE*. A file ending
  with *.og* is recommended.

Graph Sorting
-------------

| **-s, --sort**
| Apply a general topological sort to the graph and order the node ids
  accordingly. A bidirected adaptation of Kahn’s topological sort (1962)
  is used, which can handle components with no heads or tails. Here,
  both heads and tails are taken into account.

Graph Optimization
------------------
| **-O, --optimize**
| Use the MutableHandleGraph::optimize method to compact the node
  identifier space.

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for parallel operations.

Processing Information
----------------------

| **-d, --debug**
| Verbosely print graph information to stderr. This includes the maximum
  node_id, the minimum node_id, the handle to node_id mapping, the
  deleted nodes and the path metadata.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi build**.

| **-P, --progress**
| Write the current progress to stderr.

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
