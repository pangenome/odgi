.. _odgi break:

#########
odgi break
#########

break cycles in the graph and drop its paths

SYNOPSIS
========

**odgi break** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi break(1) command finds cycles in a graph via `breadth-first
search (BFS) <https://en.wikipedia.org/wiki/Breadth-first_search>`__ and
breaks them, also dropping the graph’s paths.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to break. The file name
  usually ends with *.og*.

| **-o, –out**\ =\ *FILE*
| Write the broken graph to *FILE*.

Cycle Options
-------------

| **-c, –cycle-max-bp**\ =\ *N*
| The maximum cycle length at which to break.

| **-s, –max-search-bp**\ =\ *N*
| The maximum search space of each BFS given in number of base pairs.

| **-u, –repeat-up-to**\ =\ *N*
| Iterate cycle breaking up to *N* times or stop if no new edges are
  removed.

| **-d, –show**
| Print the edges we would remove to stdout.

Program Information
-------------------

| **-h, –help**
| Print a help message for **odgi break**.

..
	EXIT STATUS
	===========

	| **0**
	| Success.

	| **1**
	| Failure (syntax or usage error; parameter error; file processing
		failure; unexpected error).
..
	BUGS
	====

	Refer to the **odgi** issue tracker at
	https://github.com/pangenome/odgi/issues.
