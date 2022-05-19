.. _odgi overlap:

#########
odgi overlap
#########

Find the paths touched by given input paths.

SYNOPSIS
========

**odgi overlap** [**-i, --input**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi overlap command finds the paths touched by the input paths.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --input**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

Overlap Options
---------------

| **-s, --subset-paths**\ =\ *FILE*
| Perform the search considering only the paths specified in the FILE;
  the file must contain one path name per line and a subset of all paths
  can be specified. When searching the overlaps, only these paths will be considered.

| **-r, --path**\ =\ *STRING*
| Perform the search of the given path *STRING* in the graph.

| **-R, --paths**\ =\ *FILE*
| Report the search results only for the paths listed in *FILE*.

| **-b, --bed-input**\ =\ *FILE*
| A BED file of ranges in paths in the graph.

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
| Print a help message for **odgi overlap**.

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
