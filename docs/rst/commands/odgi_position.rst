.. _odgi position:

#########
odgi position
#########

position parts of the graph as defined by query criteria

SYNOPSIS
========

**odgi position** [**-i, –target**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi position(1) command translates positions and coordinate ranges
between nodes and embedded paths. It provides liftover functionality,
allowing us to translate a position between any reference paths embedded
in the ``-i, --target`` graph. We can additionally project coordinates
and annotations from a source graph ``-x, --source`` into the
``target``. When completing this “graph lift”, the intersecting set of
paths in the two graphs are used to complete the coordinate projection.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –target**\ =\ *FILE*
| Describe positions in this graph.

| **-x, –source**\ =\ *FILE*
| Translate positions from this graph into the target graph using common
  ``--lift-paths`` shared between both graphs [default: use the same
  source/target graph]

Position Options
----------------

| **-r, –ref-paths**\ =\ *STRING*
| Translate the given positions into positions relative to this
  reference path.

| **-R, –ref-paths**\ =\ *FILE*
| Use the ref-paths in FILE.

| **-l, –lift-path**\ =\ *STRING*
| Lift positions from ``--source`` to ``--target`` via coordinates in
  this path common to both graphs [default: all common paths between
  ``--source`` and ``--target``].

| **-g, –graph-pos**\ =\ *[[node_id][,offset[,(+|-)]\ *\ **]**\ *]*
| A graph position, e.g. 42,10,+ or 302,0,-.

| **-F, –path-pos-file**\ =\ *FILE*
| A file with one path position per line.

| **-b, –bed-input**\ =\ *FILE*
| A BED file of ranges in paths in the graph to lift into the target
  graph ``-v``, ``--give-graph-pos`` emit graph positions.

| **-v, –give-graph-pos**
| Emit graph positions (node,offset,strand) rather than path positions.

| **-I, –all-immediate**
| Emit all positions immediately at the given graph/path position.

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
| Print a help message for **odgi position**.

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
