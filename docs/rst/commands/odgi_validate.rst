.. _odgi validate:

#############
odgi validate
#############

Validate a graph checking if the paths are consistent with the graph topology, and optionally with their source sequences.

SYNOPSIS
========

**odgi validate** [**-i, --input**\ =\ *FILE*] [**-r, --fasta**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi validate command checks that every pair of consecutive path
steps is connected by an edge in the graph.

With **-r, --fasta** it also checks that each path spells out its source
sequence, reporting counts of identical, divergent, missing path and
missing source.

OPTIONS
=======

MANDATORY OPTIONS
-----------------

| **-i, --input**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1 or GFAz (compressed GFA), but the on-the-fly conversion to the ODGI format requires additional time!

Sequence Validation
-------------------

| **-r, --fasta**\ =\ *FILE*
| Check that every path spells out its source sequence in this uncompressed FASTA *FILE*. Records pair to paths by exact name, taken up to the first whitespace. The comparison is case-insensitive but symbol-exact, so **N** does not match **A**. The first difference is reported with its source position.

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for parallel operations.

Processing Information
----------------------

| **-P, --progress**
| Print information about the operations and the progress to stderr.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi validate**.

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
