.. _odgi view:

#########
odgi view
#########

projection of graphs into other formats

SYNOPSIS
========

**odgi view** [**-i, –idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi view(1) command can convert a graph in odgi format to GFAv1. It
can reveal a graph’s internal structures for e.g. debugging processes.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to convert from. The file
  name usually ends with *.og*.

| **-g, –to-gfa**
| Write the graph in GFAv1 format to standard output.

Summary Options
---------------

| **-d, –display**
| Show the internal structures of a graph. Print to stdout the maximum
  node identifier, the minimum node identifier, the nodes vector, the
  delete nodes bit vector and the path metadata, each in a separate
  line.

Program Information
-------------------

| **-h, –help**
| Print a help message for **odgi view**.

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
