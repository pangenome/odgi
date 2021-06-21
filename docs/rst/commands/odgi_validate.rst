.. _odgi validate:

#########
odgi validate
#########

Validate a graph checking if the paths are consistent with the graph topology.

SYNOPSIS
========

**odgi validate** [**-i, --input**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi validate command validates the graph checking
if the paths are consistent with the graph topology.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --input**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

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
| Print a help message for **odgi validate**.

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
