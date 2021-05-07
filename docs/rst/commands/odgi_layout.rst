.. _odgi layout:

#########
odgi layout
#########

use SGD to make 2D layouts of the graph

SYNOPSIS
========

**odgi layout** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi layout(1) command computes 2D layouts of the graph using
stochastic gradient descent (SGD). The input graph must be sorted and
id-compacted. The algorithm itself is described in `Graph Drawing by
Stochastic Gradient Descent <https://arxiv.org/abs/1710.04626>`__. The
force-directed graph drawing algorithm minimizes the graph’s energy
function or stress level. It applies SGD to move a single pair of nodes
at a time.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to layout. The file name
  usually ends with *.og*.

| **-o, –out**\ =\ *FILE*
| Write the rendered layout in SVG format to *FILE*.

SGD Options
-----------

| **-m, –iter-max**\ =\ *N*
| The maximum number of iterations to run the layout. Default is *30*.

| **-p, –n-pivots**\ =\ *N*
| The number of pivots for sparse layout. Default is *0* leading to a
  non-sparse layout.

| **-e, –eps**\ =\ *N*
| The learning rate for SGD layout. Default is *0.01*.

SVG Options
-----------

| **-x, –x-padding**\ =\ *N*
| The padding between the connected component layouts. Default is
  *10.0*.

| **-R, –render-scale**\ =\ *N*
| SVG scaling Default is *5.0*.

Processing Information
----------------------

| **-d, –debug**
| Print information about the components to stdout.

Program Information
-------------------

| **-h, –help**
| Print a help message for **odgi layout**.

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
