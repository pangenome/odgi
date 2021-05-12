.. _odgi unitig:

#########
odgi unitig
#########

output unitigs of the graph

SYNOPSIS
========

**odgi unitig** [**-i, –idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi unitig(1) command can print all
`unitigs <https://github.com/mcveanlab/mccortex/wiki/unitig>`__ of a
given odgi graph to standard output in FASTA format. Unitigs can also be
emitted in a fixed sequence quality FASTQ format. Various parameters can
refine the unitigs to print.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to convert from. The file
  name usually ends with *.og*.

FASTQ Options
-------------

| **-f, –fake-fastq**
| Write the unitigs in FASTQ format to stdout with a fixed quality value
  of *I*.

Unitig Options
--------------

| **-t, –sample-to**\ =\ *N*
| Continue unitigs with a random walk in the graph so that they have at
  least the given *N* length.

| **-p, –sample-plus**\ =\ *N*
| Continue unitigs with a random walk in the graph by *N* past their
  natural end.

| **-l, –min-begin-node-length**\ =\ *N*
| Only begin unitigs collection from nodes which have at least length
  *N*.

Program Information
-------------------

| **-h, –help**
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
