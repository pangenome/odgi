.. _odgi validate:

#########
odgi validate
#########

validate the graph (currently, it checks if the paths
are consistent with the graph topology)

SYNOPSIS
========

**odgi validate** [**-i, –input**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi validate(1) command validates the graph (currently, it checks
if the paths are consistent with the graph topology).

OPTIONS
=======

Graph Files IO
--------------

| **-i, –input**\ =\ *FILE*
| Validate this graph.

Threading
---------

| **-t, –threads**\ =\ *N*
| Number of threads to use.

Program Information
-------------------

| **-h, –help**
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
