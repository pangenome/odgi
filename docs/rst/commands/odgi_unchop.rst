.. _odgi unchop:

#########
odgi unchop
#########

Merge unitigs into a single node preserving the node order.

SYNOPSIS
========

**odgi unchop** [**-i, --idx**\ =\ *FILE*] [**-o, --out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi unchop command merges each unitig into a single node
preserving the node order.

OPTIONS
=======

Graph Files IO
--------------

| **-i, --idx**\ =\ *FILE*
| File containing the succinct variation graph to unchop. The file name
  usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

| **-o, --out**\ =\ *FILE*
| Write the unchopped dynamic succinct variation graph in ODGI format to this *FILE*. A file ending with *.og* is recommended.

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for the parallel operations.

Processing Information
----------------------

| **-P, --progress**
| Print information about the process to stderr.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi unchop**.

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
