.. _odgi matrix:

#########
odgi matrix
#########

write the graph topology in sparse matrix formats

SYNOPSIS
========

**odgi matrix** [**-i, –idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi matrix(1) command generates a sparse matrix format out of the
graph topology of a given variation graph.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to create the sparse
  matrix from. The file name usually ends with *.og*.

Matrix Options
--------------

| **-e, –edge-depth-weight**
| Weigh edges by their path depth.

| **-d, –delta-weight**
| Weigh edges by their inverse id delta.

Program Information
-------------------

| **-h, –help**
| Print a help message for **odgi matrix**.

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
