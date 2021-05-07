.. _odgi build:

#########
odgi build
#########

construct a dynamic succinct variation graph

SYNOPSIS
========

**odgi build** [**-g, –gfa**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi build(1) command constructs a succinct variation graph from a
GFA. Currently, only GFA1 is supported. For details of the format please
see https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md.

OPTIONS
=======

Graph Files IO
--------------

| **-g, –gfa**\ =\ *FILE*
| GFA1 file containing the nodes, edges and paths to build a dynamic
  succinct variation graph from.

| **-o, –out**\ =\ *FILE*
| Write the dynamic succinct variation graph to this file. A file ending
  with *.og* is recommended.

Graph Sorting
-------------

| **-s, –sort**
| Apply a general topological sort to the graph and order the node ids
  accordingly. A bidirected adaptation of Kahn’s topological sort (1962)
  is used, which can handle components with no heads or tails. Here,
  both heads and tails are taken into account.

Processing Information
----------------------

| **-p, –progress**
| Print progress updates to stdout.

| **-d, –debug**
| Verbosely print graph information to stderr. This includes the maximum
  node_id, the minimum node_id, the handle to node_id mapping, the
  deleted nodes and the path metadata.

| **–trace**
| Include backtrace information when reporting errors.

| **-v, –verbose**
| Verbosely print processing information to stderr, including
  debug-level log messages.

| **-w, –warnings**
| Turn on script warnings (applies to executed code).

| **-t, –threads**\ =\ *N*
| Number of threads to use for the parallel operations.

Program Information
-------------------

| **-h, –help**
| Print a help message for **odgi build**.

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
