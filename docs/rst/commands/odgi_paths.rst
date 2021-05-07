.. _odgi paths:

#########
odgi paths
#########

embedded path interrogation

SYNOPSIS
========

**odgi paths** [**-i, –idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi paths(1) command allows the investigation of paths of a given
variation graph. It can calculate overlap statistics of groupings of
paths.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to investigate the paths
  from. The file name usually ends with *.og*.

| **-O, –overlaps**\ =\ *FILE*
| Read in the path grouping file to generate the overlap statistics
  from. The file must be tab-delimited. The first column lists a
  grouping and the second the path itself. Each line has one path entry.
  For each group the pairwise overlap statistics for each pairing will
  be calculated and printed to stdout.

Investigation Options
---------------------

| **-L, –list-paths**
| Print the paths in the graph to stdout. Each path is printed in its
  own line.

| **-H, –haplotypes**
| Print to stdout the paths in an approximate binary haplotype matrix
  based on the graph’s sort order. The output is tab-delimited:
  **path.name**, **path.length**, **node.count**, **node.1**,
  **node.2**, **node.n**. Each path entry is printed in its own line.

| **-D, –delim**\ =\ *CHAR*
| The part of each path name before this delimiter is a group
  identifier. This parameter should only be set in combination with
  [**-H, –haplotypes**]. Prints an additional, first column
  **group.name** to stdout.

| **-d, –distance**
| Provides a sparse distance matrix for paths. If [**-D, –delim**] is
  set, it will be path groups distances.

| **-f, –fasta**
| Print paths in FASTA format to stdout.

Threading
---------

| **-t, –threads**\ =\ *N*
| Number of threads to use.

Program Information
-------------------

| **-h, –help**
| Print a help message for **odgi paths**.

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
