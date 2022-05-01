.. _odgi inject:

#########
odgi inject
#########

Inject BED interval ranges as paths in the graph.

SYNOPSIS
========

**odgi inject** [**-i, --input**\ =\ *FILE*] [**-b, --bed-targets**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi inject command converts BED records against graph paths into new paths labeled by the BED record name.
Injection allows us to import genome annotations as paths, and is useful to produce input to odgi untangle.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --input**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

| **-o, --out**\ =\ *FILE*
| Write the graph with injected paths to *.og*.

Injection Options
----------------

| **-b, --bed-targets**\ =\ *FILE*
| BED file over path space of the graph. Records will be converted into new paths in the output graph.

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
| Print a help message for **odgi inject**.

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
