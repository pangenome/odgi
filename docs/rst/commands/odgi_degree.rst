.. _odgi degree:

#########
odgi degree
#########

describe the graph in terms of node degree

SYNOPSIS
========

**odgi degree** [**-i, –idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi degree(1) command describes the graph in terms of node degree.
For the input graph, it shows the node.count, edge.count, avg.degree,
min.degree, and max.degree.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| Describe node degree in this graph. The file name usually ends with
  *.og*.

Summary Options
---------------

| **-S, –summarize**
| Summarize the graph properties and dimensions. Print to stdout the
  node.id and the node.degree.

Threading
---------

| **-t, –threads**\ =\ *N*
| Number of threads to use.

Program Information
-------------------

| **-h, –help**
| Print a help message for **odgi degree**.

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
