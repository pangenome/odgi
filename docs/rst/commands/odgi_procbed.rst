.. _odgi procbed:

#########
odgi procbed
#########

ProcBED, or "procrustes-BED", cuts and adjusts BED intervals to fit inside a subgraph.
Coordinates of subgraph paths are taken from their names using `PanSN sequence naming format <https://github.com/pangenome/PanSN-spec>`_.
Graphs with this naming format are extracted with **odgi extract**.
**odgi procbed** can be used to map BED annotations against full reference paths to BED records against the sub-ranges of reference paths in a subgraph.
This is useful to produce input to **odgi inject**.

SYNOPSIS
========

**odgi procbed** [**-i, --input**\ =\ *FILE*] [**-b, --bed-targets**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi procbed command converts BED records against graph paths into new paths labeled by the BED record name.
Procbedion allows us to import genome annotations as paths, and is useful to produce input to odgi untangle.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --input**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

| **-o, --out**\ =\ *FILE*
| Write the graph with procbeded paths to *.og*.

Procbedion Options
----------------

| **-b, --bed-targets**\ =\ *FILE*
| BED file over path space of the full graph from the input subgraph was generated. Where they fully overlap, records will be adjusted to fit in the subgraph coordinate space.

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
| Print a help message for **odgi procbed**.

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
