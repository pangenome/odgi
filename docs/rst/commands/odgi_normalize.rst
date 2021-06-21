.. _odgi normalize:

#########
odgi normalize
#########

Compact unitigs and simplify redundant furcations.

SYNOPSIS
========

**odgi normalize** [**-i, --idx**\ =\ *FILE*] [**-o, --out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi normalize command unchops
:ref:`odgi unchop` a given variation graph
and simplifies redundant furcations.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*.

| **-o, --out**\ =\ *FILE*
| Write the normalized dynamic succinct variation graph in ODGI format to this file. A
  file ending with *.og* is recommended.

Normalize Options
-----------------

| **-I, --max-iterations**\ =\ *N*
| Iterate the normalization up to *N* many times. The default is *10*.

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for parallel operations.

Processing Information
-----------------

| **-d, --debug**
| Print information about the normalization process to stdout.

| **-P, --progress**
| Print information about the operations and the progress to stderr.

Program Information
-------------------

| **-h, --help**
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
