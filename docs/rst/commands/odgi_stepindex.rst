.. _odgi stepindex:

#########
odgi stepindex
#########

Generate a step index from a given graph. If no output file is provided via **-o, --out**, the index will be directly written to **INPUT_GRAPH.stpidx**.

SYNOPSIS
========

**odgi stepindex** [**-i, --idx**\ =\ *FILE*] [**-o, --out**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi stepindex command generates a step index from a given graph. Such an index allows us to efficiently retrieve the nucleotide position of a given graph step.
In order to save memory, a sampled step index is implemented here. We solve memory issues by only indexing every node with node identifier fitting **mod(node_id, step-index-sample-rate) == 0** in the graph.
From a given step, we can find its position by walking backwards until a node fitting our sampling criteria is found. We can retrieve this position easily, adding up the walked distance to retrieve the actual position of the step.
Effectively, the sample rate is only allowed to be a number by the power of 2, because we can use bit shift operations to calculate the modulo in O(1)! (`https://www.geeksforgeeks.org/compute-modulus-division-by-a-power-of-2-number/ <https://www.geeksforgeeks.org/compute-modulus-division-by-a-power-of-2-number/>`_).
As `evaluated <https://docs.google.com/presentation/d/1a8bOnulta6fYnQ2DFmdzt4es2vaRGmgIxO3kCe-HXR8/edit#slide=id.p>`_, the default sample rate is 8, which represents a good compromise between performance and memory usage. For ultra large graphs with hundreds of gigabytes in size, a sample rate of 16 might suite better.

As a bonus, the step index includes all the lengths of the paths, too. This allows us to efficiently get the length in nucleotides of a path by a given path handle.

Current ODGI tools that work with a step index are :ref:`odgi untangle` and :ref:`odgi tips`.

OPTIONS
=======

MANDATORY OPTIONS
-----------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

| **-o, --out**\ =\ *FILE*
| Write the created step index to the specified file. A file ending with *.stpidx* is recommended. (default: *INPUT_GRAPH.stpidx*).

Step Index Options
-------------

| **-a, --step-index-sample-rate**\ =\ *N*
| The sample rate when building the step index. We index a node only if **mod(node_id, step-index-sample-rate) == 0**! Number must be dividable by 2 or 0 to disable sampling. (default: 8).

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
| Print a help message for **odgi tips**.

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
