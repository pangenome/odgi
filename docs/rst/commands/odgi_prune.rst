.. _odgi prune:

#########
odgi prune
#########

remove complex parts of the graph

SYNOPSIS
========

**odgi prune** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi prune(1) command can remove complex parts of a graph. One can
drop paths, nodes by a certain kind of edge coverage, edges and graph
tips. Specifying a kmer length and a maximum number of furcations, the
graph can be broken at edges not fitting into these conditions.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to load in. The file name
  usually ends with *.og*.

| **-o, –out**\ =\ *FILE*
| Write the pruned graph to *FILE*. The file name should end with *.og*.

Kmer Options
------------

| **-k, –kmer-length**\ =\ *N*
| The length of the kmers to consider.

| **-e, –max-furcations**\ =\ *N*
| Break at edges that would induce *N* many furcations in a kmer.

Node Options
------------

| **-d, –max-degree**\ =\ *N*
| Remove nodes that have a higher node degree than *N*.

| **-c, –min-coverage**\ =\ *N*
| Remove nodese covered by fewer than *N* number of path steps.

| **-C, –max-coverage**\ =\ *N*
| Remove nodes covered by more than *N* number of path steps.

| **-T, –cut-tips**\ =\ *N*
| Remove nodes which are graph tips.

Edge Options
------------

| **-E, –edge-coverage**
| Remove edges outside of the minimum and maximum coverage rather than
  nodes. Only set this argument in combination with [**-c,
  –min-coverage**\ =\ *N*] and [**-C, –max-coverage**\ =\ *N*].

| **-b, –best-edges**\ =\ *N*
| Only keep the *N* most covered inbound and output edges of each node.

Path Options
------------

| **-D, –drop-paths**
| Remove the paths from the graph.

Threading
---------

| **-t, –threads**\ =\ *N*
| Number of threads to use.

Program Information
-------------------

| **-h, –help**
| Print a help message for **odgi prune**.

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
