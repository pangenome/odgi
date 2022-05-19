.. _odgi flatten:

#########
odgi flatten
#########

Generate linearizations of a graph.

SYNOPSIS
========

**odgi flatten** [**-i, --idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi flatten command projects the graph sequence and paths into
FASTA and BED.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

Output Options
--------------

| **-f, --fasta**\ =\ *FILE*
| Write the concatenated node sequences in FASTA format to *FILE*.

| **-n, --name-seq**\ =\ *STRING*
| The name to use for the concatenated graph sequence. Default is the
  name of the input file which was specified via [**-i,
  --idx**\ =\ *FILE*].

| **-b, --bed**\ =\ *FILE*
| Write the mapping between graph paths and the linearized FASTA
  sequence in BED format to *FILE*.

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use (to embed the subpaths in parallel).

Processing Information
----------------------

| **-P, --progress**
| Print information about the operations and the progress to stderr.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi flatten**.

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
