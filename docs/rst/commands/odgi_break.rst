.. _odgi break:

#########
odgi break
#########

Break cycles in the graph and drop its paths.

SYNOPSIS
========

**odgi break** [**-i, --idx**\ =\ *FILE*] [**-o, --out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi break command finds cycles in a graph via `breadth-first
search (BFS) <https://en.wikipedia.org/wiki/Breadth-first_search>`__ and
breaks them, also dropping the graph’s paths.

OPTIONS
=======

MANDATORY OPTIONS
-----------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*.

| **-o, --out**\ =\ *FILE*
| Write the broken graph in ODGI format to *FILE*. A file ending of *.og* is recommended.

Cycle Options
-------------

| **-c, --cycle-max-bp**\ =\ *N*
| The maximum cycle length at which to break (default: 0).

| **-s, --max-search-bp**\ =\ *N*
| The maximum search space of each BFS given in number of base pairs (default: 0).

| **-u, --repeat-up-to**\ =\ *N*
| Iterate cycle breaking up to *N* times or stop if no new edges are
  removed.

| **-d, --show**
| Print the edges we would remove to stdout.

Program Information
-------------------

| **-h, --help**
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
