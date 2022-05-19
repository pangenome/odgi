.. _odgi explode:

#########
odgi explode
#########

Breaks a graph into connected components storing each component in its own file.

SYNOPSIS
========

**odgi explode** [**-i, --idx**\ =\ *FILE*] [**-p,
–prefix**\ =\ *STRING*] [*OPTION*]…

DESCRIPTION
===========

The odgi explode command breaks a graph into connected components,
writing each component in its own file.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

Explode Options
---------------

| **-g, --to-gfa**
| Write each connected component to a file in GFAv1 format.

| **-p, --prefix**\ =\ *STRING*
| Write each connected component in a file with the given *STRING* prefix. The
  file for the component number ``i`` will be named ``STRING.i.EXTENSION``
(default: ``component.i.og`` or ``component.i.gfa``).

| **-b, --biggest**\ =\ *N*
| Specify the number of the biggest connected components to write,
  sorted by decreasing size (default: disabled, for writing them all).

| **-s, --sorting-criteria**\ =\ *C*
| Specify how to sort the connected components by size:
| p) Path mass (total number of path bases) (default).
| l) Graph length (number of node bases).
| n) Number of nodes.
| P) Longest path.

| **-O, --optimize**
| Compact the node ID space in each connected component.

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for parallel operations.

Processing Information
----------------------

| **-P, --progress**
| Print information about the components and the progress to stderr.

Program Information
-------------------

| **-h, --help**
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
