.. _odgi squeeze:

#########
odgi squeeze
#########

Squeezes multiple graphs in ODGI format into the same file in ODGI format.

SYNOPSIS
========

**odgi squeeze** [**-f, --input-graphs**\ =\ *FILE*] [**-o,
–out**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi squeeze command merges multiple graphs into the same file.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-f, --input-graphs**\ =\ *FILE*
| Input file containing the list of graphs to squeeze into the same
  file. The file must contain one graph per line. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

| **-o, --out**\ =\ *FILE*
| Store all the input graphs in this *FILE*. The file name usually ends with *.og*.

Squeeze Options
---------------

| **-s, --rank-suffix**\ =\ *STRING*
| Add the separator and the input file rank as suffix to the path names
  (to avoid path name collisions).

| **-O, --optimize**
| Compact the node ID space for each connected component before squeezing.

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for parallel operations.

Processing Information
----------------------

| **-P, --progress**
| Print information about the progress to stderr.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi squeeze**.

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
