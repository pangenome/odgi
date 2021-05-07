.. _odgi explode:

#########
odgi explode
#########

breaks a graph into connected components in their own
files

SYNOPSIS
========

**odgi explode** [**-i, –idx**\ =\ *FILE*] [**-p,
–prefix**\ =\ *STRING*] [*OPTION*]…

DESCRIPTION
===========

The odgi explode(1) command breaks a graph into connected components,
writing each component in its own file.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to break in its
  components. The file name usually ends with *.og*.

Explode Options
---------------

| **-p, –prefix**\ =\ *STRING*
| Write each connected component in a file with the given prefix. The
  file for the component ``i`` will be named ``STRING.i.og`` (default:
  ``component``).

| **-b, –biggest**\ =\ *N*
| Specify the number of the biggest connected components to write,
  sorted by decreasing size (default: disabled, for writing them all).

| **-s, –sorting-criteria**\ =\ *C*
| Specify how to sort the connected components by size:

-  p) path mass (total number of path bases) (default)

-  l) graph length (number of node bases)

-  n) number of nodes

-  P) longest path

| **-O, –optimize**
| Compact the node ID space in each connected component.

Threading
---------

| **-t, –threads**\ =\ *N*
| Number of threads to use (to write the components in parallel).

Processing Information
----------------------

| **-P, –progress**
| Print information about the components and the progress to stderr.

Program Information
-------------------

| **-h, –help**
| Print a help message for **odgi explode**.

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
