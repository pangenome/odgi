.. _odgi unchop:

#########
odgi unchop
#########

merge unitigs into single nodes

SYNOPSIS
========

**odgi unchop** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi unchop(1) command merges each unitig into a single node
preserving the node order.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to unchop. The file name
  usually ends with *.og*.

| **-o, –out**\ =\ *FILE*
| Write the unchopped dynamic succinct variation graph to this file. A
  file ending with *.og* is recommended.

Processing Information
----------------------

| **-d, –debug**
| Print information about the process to stderr.

Threading
---------

| **-t, –threads**\ =\ *N*
| Number of threads to use for the parallel operations.

Program Information
-------------------

| **-h, –help**
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
