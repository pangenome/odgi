.. _odgi panpos:

#########
odgi panpos
#########

Get the pangenome position of a given path and nucleotide position (1-based).

SYNOPSIS
========

**odgi panpos** [**-i, --idx**\ =\ *FILE*] [**-p, --path**\ =\ *STRING*]
[**-n, --nuc-pos**\ =\ *N*] [*OPTION*]…

DESCRIPTION
===========

The odgi panpos command give a pangenome position for a given path
and nucleotide position. It requires a path index, which can be created
with :ref:`odgi pathindex`. Going
from **path:position** → **pangenome:position** is important when
navigating large graphs in an interactive manner like in the
`Pantograph <https://graph-genome.github.io/>`__ project. All input and
output positions are 1-based.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph index in xp format from this *FILE*. The file name usually ends with *.xp*.

| **-p, --path**\ =\ *STRING*
| The path name of the query.

| **-n, --nuc-pos**\ =\ *STRING*
| The nucleotide query of the query.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi panpos**.

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
