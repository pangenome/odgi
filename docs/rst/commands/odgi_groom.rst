.. _odgi groom:

#########
odgi groom
#########

Resolve spurious inverting links.

SYNOPSIS
========

**odgi groom** [**-i, --idx**\ =\ *FILE*] [**-o, --out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi groom command resolves spurious inverting links.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*.

| **-o, --out**\ =\ *FILE*
| Write the groomed succinct variation graph in ODGI format to *FILE*. A file ending with *.og* is recommended.

Grooming Options
----------------

| **-d, --use-dfs**
| Use depth-first search for the grooming.

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for parallel operations.

Processing Information
----------------------

| **-P, --progress**
| Write the current progress to stderr.

Program Information
-------------------

| **-h, --help**
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
