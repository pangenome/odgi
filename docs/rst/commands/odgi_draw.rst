.. _odgi draw:

#########
odgi draw
#########

Draw previously-determined 2D layouts of the graph with diverse annotations.

SYNOPSIS
========

**odgi draw** [**-i, --idx**\ =\ *FILE*] [**-c, --coords-in**\ =\ *FILE*]
[**-p, --png**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi draw command draws previously-determined 2D layouts of the
graph with diverse annotations.

OPTIONS
=======

MANDATORY OPTIONS
-----------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

| **-c, --coords-in**\ =\ *FILE*
| Read the layout coordinates from this .lay format *FILE* produced by :ref:`odgi layout`.

Files IO
--------

| **-T, --tsv**\ =\ *FILE*
| Write the TSV layout plus displayed annotations to this *FILE*.

| **-s, --svg**\ =\ *FILE*
| Write an SVG rendering to this *FILE*.

| **-p, --png**\ =\ *FILE*
| Write a rasterized PNG rendering to this *FILE*.

| **-X, --path-index**\ =\ *FILE*
| Load the path index from this *FILE*.

Visualization Options
---------------------

| **-H, --png-height**\ =\ *N*
| Height of PNG rendering (default: 1000).

| **-E, --png-border**\ =\ *N*
| Size of PNG border in bp (default: 10).

| **-C –color-paths**
| Color paths (in PNG output).

| **-R, --scale**\ =\ *N*
| Image scaling (default 1.0).

| **-B, --border**\ =\ *N*
| Image border (in approximate bp) (default 100.0).

| **-w, --line-width**\ =\ *N*
| Line width (in approximate bp) (default 0.0).

| **-S, --path-line-spacing**\ =\ *N*
| Spacing between path lines in PNG layout (in approximate bp) (default
  0.0).

| **-b, --bed-file**\ =\ *FILE*
Color the nodes based on the input annotation in the given BED FILE.
Colors are derived from the 4th column, if present, else from the path name.
If the 4th column value is in the format 'string#RRGGBB', the RRGGBB color (in hex notation) will be used.

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
| Print a help message for **odgi draw**.

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
	
	::
	
	   Refer to the *odgi* issue tracker at https://github.com/pangenome/odgi/issues.
