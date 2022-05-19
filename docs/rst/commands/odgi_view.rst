.. _odgi view:

#########
odgi view
#########

Project a graph into other formats.

SYNOPSIS
========

**odgi view** [**-i, --idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi view command can convert a graph in odgi format to GFAv1. It
can reveal a graph’s internal structures for e.g. debugging processes.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

Output Options
--------------

| **-g, --to-gfa**
| Write the graph in GFAv1 format to standard output.

| **-a, --node-annotation**
| Emit node annotations for the graph in GFAv1 format.

Summary Options
---------------

| **-d, --display**
| Show the internal structures of a graph. Print to stderr the maximum
  node identifier, the minimum node identifier, the nodes vector, the
  delete nodes bit vector and the path metadata, each in a separate
  line.

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
