.. _odgi draw:

#########
odgi draw
#########

variation graph visualizations in 2D

SYNOPSIS
========

**odgi draw** [**-i, –idx**\ =\ *FILE*] [**-c, –coords-in**\ =\ *FILE*]
[**-p, –png**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi draw(1) command draws previously-determined 2D layouts of the
graph with diverse annotations.

OPTIONS
=======

Files IO
--------

| **-i, –idx**\ =\ *FILE*
| Load the graph from this file. The file name usually ends with *.og*.

| **-c, –coords-in**\ =\ *FILE*
| Read the layout coordinates from this .lay format file produced by :ref:`odgi layout`.

| **-T, –tsv**\ =\ *FILE*
| Write the TSV layout plus displayed annotations to this file.

| **-s, –svg**\ =\ *FILE*
| Write an SVG rendering to this file.

| **-p, –png**\ =\ *FILE*
| Write a rasterized PNG rendering to this file.

Visualization Options
---------------------

| **-H, –png-height**\ =\ *N*
| Height of PNG rendering (default: 1000).

| **-E, –png-border**\ =\ *N*
| Size of PNG border in bp (default: 10).

| **-C –color-paths**
| Color paths (in PNG output).

| **-R, –scale**\ =\ *N*
| Image scaling (default 1.0).

| **-B, –border**\ =\ *N*
| Image border (in approximate bp) (default 100.0).

| **-w, –line-width**\ =\ *N*
| Line width (in approximate bp) (default 0.0).

| **-S, –path-line-spacing**\ =\ *N*
| Spacing between path lines in png layout (in approximate bp) (default
  0.0).

| **-X, –path-index**\ =\ *FILE*
| Load the path index from this file.

Program Information
-------------------

| **-h, –help**
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
