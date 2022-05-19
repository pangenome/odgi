.. _odgi pathindex:

#########
odgi pathindex
#########

Create a path index for a given graph.

SYNOPSIS
========

**odgi pathindex** [**-i, --idx**\ =\ *FILE*] [**-o, --out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi pathindex command generates a path index of a graph. It uses
succinct data structures to encode the index. The path index represents
a subset of the features of a fully realized `xg
index <https://github.com/vgteam/xg>`__. Having a path index, we can use
:ref:`odgi panpos` to go from
**path:position** → **pangenome:position** which is important when
navigating large graphs in an interactive manner like in the
`Pantograph <https://graph-genome.github.io/>`__ project.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

| **-o, --out**\ =\ *FILE*
| Write the succinct variation graph index to this FILE. A file ending with *.xp* is recommended.

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
| Print a help message for **odgi pathindex**.

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
