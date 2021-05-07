.. _odgi groom:

#########
odgi groom
#########

resolve spurious inverting links

SYNOPSIS
========

**odgi groom** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi groom(1) command resolves spurious inverting links.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to groom. The file name
  usually ends with *.og*.

| **-o, –out**\ =\ *FILE*
| Write the groomed succinct variation graph to *FILE*. The file name
  usually ends with *.og*.

Processing Information
----------------------

| **-P, –progress**
| Display progress of the grooming to stderr.

Program Information
-------------------

| **-h, –help**
| Print a help message for **odgi groom**.

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
