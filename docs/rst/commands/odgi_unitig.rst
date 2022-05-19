.. _odgi unitig:

#########
odgi unitig
#########

Output unitigs of the graph.

SYNOPSIS
========

**odgi unitig** [**-i, --idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi unitig command can print all
`unitigs <https://github.com/mcveanlab/mccortex/wiki/unitig>`__ of a
given odgi graph to standard output in FASTA format. Unitigs can also be
emitted in a fixed sequence quality FASTQ format. Various parameters can
refine the unitigs to print.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

FASTQ Options
-------------

| **-f, --fake-fastq**
| Write the unitigs in FASTQ format to stdout with a fixed quality value of *I*.

Unitig Options
--------------

| **-t, --sample-to**\ =\ *N*
| Continue unitigs with a random walk in the graph so that they have at least the given *N* length.

| **-p, --sample-plus**\ =\ *N*
| Continue unitigs with a random walk in the graph by *N* past their natural end.

| **-l, --min-begin-node-length**\ =\ *N*
| Only begin unitigs collection from nodes which have at least length *N*.

Processing Information
----------------------

| **-P, --progress**
| Print information about the operations and the progress to stderr.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi extract**.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi unitig**.

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
