.. _odgi similarity:

#########
odgi similarity
#########

Provides a sparse similarity matrix for paths or groups of paths.
Each line prints in a tab-delimited format to stdout.

SYNOPSIS
========

**odgi similarity** [**-i, --idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi similarity command allows the investigation of the similarity between (groups of) paths of a given variation graph.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

Path Investigation Options
---------------------

| **-D, --delim**\ =\ *CHAR*
| The part of each path name before this delimiter is a group identifier.

| **-p, --delim-pos**\ =\ *N*
| Consider the N-th occurrence of the delimiter specified with **-D, --delim** to obtain the
  group identifier. Specify 1 for the 1st occurrence (default)."

| **-d, --distances**
| Provide distances (dissimilarities) instead of similarities.
  Outputs an additional column with the Euclidean distance.

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
| Print a help message for **odgi paths**.

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
