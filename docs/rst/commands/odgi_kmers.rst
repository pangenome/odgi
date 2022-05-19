.. _odgi kmers:

#########
odgi kmers
#########

Display and characterize the kmer space of a graph.

SYNOPSIS
========

**odgi kmers** [**-i, --idx**\ =\ *FILE*] [**-c, --stdout**] [*OPTION*]…

DESCRIPTION
===========

Given a kmer length, the odgi kmers command can emit all kmers. The
output can be refined by setting the maximum number of furcations at
edges or by not considering nodes above a given node degree limit.

OPTIONS
=======

MANDATORY ARGUMENTS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

| **-k, --kmer-length**\ =\ *N*
| The kmer length to generate kmers from.

Kmer Options
------------

| **-c, --stdout**
| Write the kmers to standard output. Kmers are line-separated.

| **-e, --max-furcations**\ =\ *N*
| Break at edges that would induce this many furcations when generating
  a kmer.

| **-D, --max-degree**\ =\ *N*
| Don’t take nodes into account that have a degree greater than *N*.

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
| Print a help message for **odgi kmers**.

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
