.. _odgi pathindex:

#########
odgi pathindex
#########

create a path index for a given path

SYNOPSIS
========

**odgi pathindex** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi pathindex(1) command generates a path index of a graph. It uses
succinct data structures to encode the index. The path index represents
a subset of the features of a fully realized `xg
index <https://github.com/vgteam/xg>`__. Having a path index, we can use
:ref:`odgi panpos` to go from
**path:position** → **pangenome:position** which is important when
navigating large graphs in an interactive manner like in the
`Pantograph <https://graph-genome.github.io/>`__ project.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to generate a path index
  from. The file name usually ends with *.og*.

| **-o, –out**\ =\ *FILE*
| Write the path index to *FILE*.

Program Information
-------------------

| **-h, –help**
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
