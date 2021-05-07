.. _odgi chop:

#########
odgi chop
#########

divide nodes into smaller pieces

SYNOPSIS
========

**odgi chop** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*] [**-c,
–chop-to**\ =\ *N*] [*OPTION*]…

DESCRIPTION
===========

The odgi chop(1) command chops long nodes into short ones while
preserving the graph topology and node order.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to chop. The file name
  usually ends with *.og*.

| **-o, –out**\ =\ *FILE*
| Write the choped succinct variation graph to *FILE*. The file name
  usually ends with *.og*.

Chop Options
------------

| **-c, –chop-to**\ =\ *N*
| Divide nodes that longer than *N* into nodes no longer than *N* while
  maintaining graph topology.

Threading
---------

| **-t, –threads**\ =\ *N*
| Number of threads to use for the parallel operations.

Processing Information
----------------------

| **-d, –debug**
| Print information about the process to stderr.

Program Information
-------------------

| **-h, –help**
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
