.. _odgi overlap:

#########
odgi overlap
#########

find the paths touched by the input paths

SYNOPSIS
========

**odgi overlap** [**-i, –input**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi overlap(1) command finds the paths touched by the input paths.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –input**\ =\ *FILE*
| Perform the search in this graph.

Overlap Options
---------------

| **-s, –subset-paths**\ =\ *FILE*
| Perform the search considering only the paths specified in the FILE;
  the file must contain one path name per line and a subset of all paths
  can be specified.

| **-r, –path**\ =\ *STRING*
| Perform the search of the given path in the graph.

| **-R, –paths**\ =\ *FILE*
| Perform the search for the paths listed in FILE

| **-b, –bed-input**\ =\ *FILE*
| A BED file of ranges in paths in the graph.

Threading
---------

| **-t, –threads**\ =\ *N*
| Number of threads to use.

Program Information
-------------------

| **-h, –help**
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
