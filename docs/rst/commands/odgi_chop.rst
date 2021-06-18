.. _odgi chop:

#########
odgi chop
#########

Divide nodes into smaller pieces preserving node topology and order.

SYNOPSIS
========

**odgi chop** [**-i, --idx**\ =\ *FILE*] [**-o, --out**\ =\ *FILE*] [**-c,
–chop-to**\ =\ *N*] [*OPTION*]…

DESCRIPTION
===========

The odgi chop command chops long nodes into short ones while
preserving the graph topology and node order.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

| **-o, --out**\ =\ *FILE*
| Write the chopped succinct variation graph in ODGI format to *FILE*. A file ending of *.og* is recommended.

| **-c, --chop-to**\ =\ *N*
| Divide nodes that are longer than *N* base pairs into nodes no longer than *N* while
  maintaining graph topology.

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for the parallel operations.

Processing Information
----------------------

| **-d, --debug**
| Print information about the process to stderr.

| **-P, --progress**
| Print information about the operations and the progress to stderr.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi chop**.

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
