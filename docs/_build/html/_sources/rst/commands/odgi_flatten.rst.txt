.. _odgi flatten:

#########
odgi flatten
#########

generate linearization of the graph

SYNOPSIS
========

**odgi flatten** [**-i, –idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi flatten(1) command projects the graph sequence and paths into
FASTA and BED.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to flatten. The file name
  usually ends with *.og*.

Output Options
--------------

| **-f, –fasta**\ =\ *FILE*
| Write the concatenated node sequences in FASTA format to *FILE*.

| **-n, –name-seq**\ =\ *STRING*
| The name to use for the concatenated graph sequence. Default is the
  name of the input file which was specified via [**-i,
  –idx**\ =\ *FILE*].

| **-b, –bed**\ =\ *FILE*
| Write the mapping between graph paths and the linearized FASTA
  sequence in BED format to *FILE*.

Program Information
-------------------

| **-h, –help**
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
