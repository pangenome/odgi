.. _odgi position:

#########
odgi position
#########

Find, translate, and liftover graph and path positions between graphs. Results are printed to stdout.

SYNOPSIS
========

**odgi position** [**-i, --target**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi position command translates positions and coordinate ranges
between nodes and embedded paths. It provides liftover functionality,
allowing us to translate a position between any reference paths embedded
in the ``-i, --target`` graph. We can additionally project coordinates
and annotations from a source graph ``-x, --source`` into the
``target``. When completing this “graph lift”, the intersecting set of
paths in the two graphs are used to complete the coordinate projection.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --target**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

Position Options
----------------

| **-x, --source**\ =\ *FILE*
| Translate positions from this *FILE* graph into the target graph using common
  **-l, --lift-paths** shared between both graphs (default: use the same
  source/target graph). It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

| **-r, --ref-path**\ =\ *STRING*
| Translate the given positions into positions relative to this
  reference path *STRING*.

| **-R, --ref-paths**\ =\ *FILE*
| Use the ref-paths in *FILE* for positional translation.

| **-l, --lift-path**\ =\ *STRING*
| Lift positions from *-x, --source* to *-i, --target* via coordinates in
  this path common to both graphs (default: all common paths between
  *-x, --source* and *-i, --target*).

| **-L, --lift-paths**\ =\ *FILE*
| Same as in *-l, --lift-paths*, but for all paths in *FILE*.

| **-g, --graph-pos**\ =\ *[[node_id][,offset[,(+|-)]\ *\ **]**\ *]*
| A graph position, e.g. 42,10,+ or 302,0,-.

| **-G, --graph-pos-file**\ =\ *FILE*
| Same as in *-g, --graph-pos*, but for all graph positions in *FILE*.

| **-p, --path-pos**\ =\ *[path_name][,offset[,(+|-)]*]**
| A path position, e.g. chr8,1337,+ or chrZ,3929,-.

| **-F, --path-pos-file**\ =\ *FILE*
| A *FILE* with one path position per line.

| **-b, --bed-input**\ =\ *FILE*
| A BED file of ranges in paths in the graph to lift into the target
  graph *-v, --give-graph-pos* emit graph positions.

| **-E, --gff-input**\ =\ *FILE*
| A GFF/GTF file with annotation of ranges in paths in the graph to lift into the target (sub)graph emitting graph identifiers with annotation. The output is a CSV reading for the visualization within Bandage. The first column is the node identifier, the second column the annotation. If several annotations exist for the same node, they are combined via ';'.


| **-v, --give-graph-pos**
| Emit graph positions (node, offset, strand) rather than path positions.

| **-I, --all-immediate**
| Emit all positions immediately at the given graph/path position.

| **-d, --search-radius**\ =\ *DISTANCE*
| Limit coordinate conversion breadth-first search up to DISTANCE bp
  from each given position (default: 10000).

| **-w, --jaccard-context**\ =\ *N*
| Maximum walking distance in nucleotides for one orientation when finding the best target (reference) range for each query path (default: 10000). Note: If we walked 9999 base pairs and **w, --jaccard-context** is **10000**, we will also include the next node, even if we overflow the actual limit.

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
