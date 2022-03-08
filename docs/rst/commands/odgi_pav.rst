.. _odgi pav:

#########
odgi pav
#########

Presence/absence variants (PAVs).
It prints to stdout a matrix with the PAVs ratios.

SYNOPSIS
========

**odgi pav** [**-i, --idx**\ =\ *FILE*] [**-b,
–bed-file**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi pav command prints to stdout a matrix with the Presence/absence variants (PAVs) ratios.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

| **-b, --bed-file**\ =\ *FILE*
| Find PAVs in the path range(s) specified in the given BED FILE

Pav Options
---------------

| **-p, --path-groups**\ =\ *FILE*
| Group paths as described in two-column *FILE*, with columns `path.name` and `group.name`.

| **-B, --binary-matrix**\ =\ *THRESHOLD*
| Emit a binary matrix, with 1 if the PAV ratio is greater than or equal to the specified THRESHOLD, else 0.

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for parallel operations.

Processing Information
----------------------

| **-P, --progress**
| Print information about the components and the progress to stderr.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi pav**.

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