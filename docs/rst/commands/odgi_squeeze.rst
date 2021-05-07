.. _odgi squeeze:

#########
odgi squeeze
#########

squeezes multiple graphs into the same file

SYNOPSIS
========

**odgi squeeze** [**-f, –input-graphs**\ =\ *FILE*] [**-o,
–out**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi squeeze(1) command merges multiple graphs into the same file.

OPTIONS
=======

Graph Files IO
--------------

| **-f, –input-graphs**\ =\ *FILE*
| Input file containing the list of graphs to squeeze into the same
  file. The file must contain one path per line.

| **-o, –out**\ =\ *FILE*
| Store all the input graphs in this file. The file name usually ends
  with *.og*.

Squeeze Options
---------------

| **-s, –rank-suffix**\ =\ *STRING*
| Add the separator and the input file rank as suffix to the path names
  (to avoid path name collisions).

| **-O, –optimize**
| Compact the node ID space in each input file before imploding.

Threading
---------

| **-t, –threads**\ =\ *N*
| Number of threads to use.

Processing Information
----------------------

| **-P, –progress**
| Print information about the progress to stderr.

Program Information
-------------------

| **-h, –help**
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
