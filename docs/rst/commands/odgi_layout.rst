.. _odgi layout:

#########
odgi layout
#########

Establish 2D layouts of the graph using path-guided stochastic gradient descent (the graph must be sorted and id-compacted).

SYNOPSIS
========

**odgi layout** [**-i, --idx**\ =\ *FILE*] [**-o, --out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi layout command computes 2D layouts of the graph using path-guided
stochastic gradient descent (PG-SGD). The input graph must be sorted and
id-compacted. The algorithm itself is described in `Graph Drawing by
Stochastic Gradient Descent <https://arxiv.org/abs/1710.04626>`__. The
force-directed graph drawing algorithm minimizes the graph’s energy
function or stress level. It applies SGD to move a single pair of nodes
at a time.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

Files IO
-------

| **-o, --out**\ =\ *FILE*
| Write the layout coordinates to this *FILE* in .lay binary format.

| **-T, --tsv**\ =\ *FILE*
| Write the layout in TSV format to this *FILE*.

| **-X, --path-index**\ =\ *FILE*
| Load the path index from this FILE so that it does
| not have to be created for the layout calculation.

| **-C, --temp-dir**\ =\ *PATH*
| Directory for temporary files.

| **-f, --path-sgd-use-paths**\ =\ *FILE*
| Specify a line separated list of paths to sample from for the on the fly term generation process in the path guided 2D SGD (default: sample from all paths).

Layout Initialization Options
-----------------------------

| **-N, --layout-initialization**\ =\ *STRING*
| Specify the layout initialization mode:
| *d*) Node rank in X and gaussian noise in Y (default).
| *r*) Uniform noise in X and Y in the order of the graph length.
| *u*) Node rank in X and uniform noise in Y.
| *g*) Gaussian noise in X and Y
| *h*) Hilbert curve in X and Y.

PG-SGD Options
--------------

| **-G, --path-sgd-min-term-updates-paths**\ =\ *N*
| Minimum number of terms *N* to be updated before a new path guided 2D SGD iteration with adjusted learning rate eta starts
 , expressed as a multiple of total path length (default: 10).

| **-U, --path-sgd-min-term-updates-nodes**\ =\ *N*
| Minimum number of terms *N* to be updated before a new path guided 2D SGD iteration with adjusted learning rate
 eta starts, expressed as a multiple of the number of nodes (default: argument is not set, the default of -G=[N],
 path-sgd-min-term-updates-paths=[N] is used.

| **-j, --path-sgd-delta**\ =\ *N*
| The threshold of the maximum displacement *N* approximately in bp at which to stop path guided 2D SGD (default: 0).

| **-g, --path-sgd-eta**\ =\ *N*
| The final learning rate *N* for path guided 2D SGD model (default: 0.01).

| **-v, --path-sgd-eta-max**\ =\ *N*
| The first and maximum learning rate *N* for path guided 2D SGD model (default: squared longest path length).

| **-a, --path-sgd-zipf-theta**\ =\ *N*
| The theta value *N* for the Zipfian distribution which is used as the sampling method for the second node of one term in
 the path guided 2D SGD model (default: 0.99).

| **-x, --path-sgd-iter-max**\ =\ *N*
| The maximum number of iterations *N* for the path guided 2D SGD model (default: 30).

| **-F, --path-sgd-iteration-max-learning-rate**\ =\ *N*
| Specify the iteration *N* where the learning rate is max for path guided 2D SGD model (default: 0).

| **-k, --path-sgd-zipf-space**\ =\ *N*
| The maximum space size *N* of the Zipfian distribution which is used as the sampling method for the second node of one
 term in the path guided 2D SGD model (default: max path lengths).

| **-I, --path-sgd-zipf-space-max**\ =\ *N*
| The maximum space size *N* of the Zipfian distribution beyond which quantization occurs (default: 1000).

| **-l, --path-sgd-zipf-space-quantization-step**\ =\ *N*
| The size of the quantization step *N* when the maximum space size of the Zipfian distribution is exceeded (default: 100).

| **-u, --path-sgd-snapshot**\ =\ *STRING*
| Set the prefix *STRING* to which each snapshot layout of a path guided 2D SGD iteration should be written to (default: NONE).

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for parallel operations.

Processing Information
----------------------

| **-P, --progress**
| Write the current progress to stderr.

Program Information
-------------------

| **-h, --help**
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
