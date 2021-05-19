.. _odgi kmers:

#########
odgi kmers
#########

show and characterize the kmer space of the graph

SYNOPSIS
========

**odgi kmers** [**-i, –idx**\ =\ *FILE*] [**-c, –stdout**] [*OPTION*]…

DESCRIPTION
===========

Given a kmer length, the odgi kmers(1) command can emit all kmers. The
output can be refined by setting the maximum number of furcations at
edges or by not considering nodes above a given node degree limit.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to convert from. The file
  name usually ends with *.og*.

| **-c, –stdout**\ =
| Write the kmers to standard output. Kmers are line-separated.

Kmer Options
------------

| **-k, –kmer-length**\ =\ *N*
| The kmer length to generate kmers from.

| **-e, –max-furcations**\ =\ *N*
| Break at edges that would induce this many furcations when generating
  a kmer.

| **-D, –max-degree**\ =\ *N*
| Don’t take nodes into account that have a degree greater than *N*.

Threading
---------

| **-t, –threads**\ =\ *N*
| Number of threads to use.

Program Information
-------------------

| **-h, –help**
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
