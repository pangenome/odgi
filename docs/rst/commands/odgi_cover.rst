.. _odgi cover:

#########
odgi cover
#########

Cover the graph with paths.

SYNOPSIS
========

**odgi cover** [**-i, --idx**\ =\ *FILE*] [**-o, --out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi cover command creates a path coverage of a variation graph, with
a specified number of paths per component.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

| **-o, --out**\ =\ *FILE*
| Write the succinct variation graph with the generated paths in ODGI format to *FILE*. A file ending with *.og* is recommended.

Cover Options
-------------

| **-H, --hogwild-depth**\ =\ *DEPTH*
| Randomly cover the graph until it reaches the given average *DEPTH*. Specifying
  this options ignores all other Cover Options except **-I, --ignore-paths**!

| **-n, --num-paths-per-component**\ =\ *N*
| Number of paths to generate per component.

| **-k, --node-windows-size**\ =\ *N*
| Size of the node window to check each time a new path is extended (it
  has to be greater than or equal to 2).

| **-c, --min-node-coverage**\ =\ *N*
| Minimum node coverage to reach (it has to be greater than 0). There
  will be generated paths until the specified minimum node coverage is
  reached, or until the maximum number of allowed generated paths is
  reached (number of nodes in the input variation graph).

| **-I, --ignore-paths**
| Ignore the paths already embedded in the graph during the nodes
  coverage initialization.

| **-w, --write-node-coverages**\ =\ *FILE*
| Write the node coverages at the end of the paths generation to *FILE*.
  The file will contain tab-separated values (header included), and have
  3 columns: *component_id*, *node_id*, and *coverage*.

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
| Print a help message for **odgi cover**.

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
