.. _odgi normalize:

#########
odgi normalize
#########

compact unitigs and simplify redundant furcations

SYNOPSIS
========

**odgi normalize** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi normalize(1) command unchops
:ref:`odgi unchop` a given variation graph
and simplifies redundant furcations.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to normalize. The file
  name usually ends with *.og*.

| **-o, –out**\ =\ *FILE*
| Write the normalized dynamic succinct variation graph to this file. A
  file ending with *.og* is recommended.

| **-I, –max-iterations**\ =\ *N*
| Iterate the normalization up to *N* many times. The default is *10*.

Program Debugging
-----------------

| **-d, –debug**
| Print information about the normalization process to stdout.

Program Information
-------------------

| **-h, –help**
| Print a help message for **odgi normalize**.

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
