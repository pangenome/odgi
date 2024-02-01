.. _odgi sort:

#########
odgi sort
#########

Apply different kinds of sorting algorithms to a graph. The most prominent one is the PG-SGD sorting algorithm.

SYNOPSIS
========

**odgi sort** [**-i, --idx**\ =\ *FILE*] [**-o, --out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi sort command sorts a succinct variation graph. Odgi sort
offers a diverse palette of sorting algorithms to determine the node
order:

-  A topological sort: A graph can be sorted via `breadth-first search
   (BFS) <https://en.wikipedia.org/wiki/Breadth-first_search>`__ or
   `depth-first search
   (DFS) <https://en.wikipedia.org/wiki/Depth-first_search>`__.
   Optionally, a chunk size specifies how much of the graph to grab at
   once in each topological sorting phase. The sorting algorithm will
   continue the sort from the next node in the prior graph order that
   has not been sorted, yet. The cycle breaking algorithm applies a DFS
   sort until a cycle is found. We break and start a new DFS sort phase
   from where we stopped.

-  A random sort: The graph is randomly sorted. The node order is
   randomly shuffled from `Mersenne Twister
   pseudo-random <http://www.cplusplus.com/reference/random/mt19937/>`__
   generated numbers.

-  A 1D linear SGD sort: ODGI implements a 1D linear, variation graph
   adjusted, multi-threaded version of the `Graph Drawing by Stochastic
   Gradient Descent <https://arxiv.org/abs/1710.04626>`__ algorithm. The
   force-directed graph drawing algorithm minimizes the graph’s energy
   function or stress level. It applies stochastic gradient descent
   (SGD) to move a single pair of nodes at a time.

-  A path guided, 1D linear SGD sort: ODGI implements a 1D linear,
   variation graph adjusted, multi-threaded version of the `Graph
   Drawing by Stochastic Gradient
   Descent <https://arxiv.org/abs/1710.04626>`__ algorithm. The
   force-directed graph drawing algorithm minimizes the graph’s energy
   function or stress level. It applies stochastic gradient descent
   (SGD) to move a single pair of nodes at a time. The path index is
   used to pick the terms to move stochastically. For more details about
   the algorithm, please take a look at https://www.biorxiv.org/content/10.1101/2023.09.22.558964v2.

Sorting the paths in a graph my refine the sorting process. For the
users’ convenience, it is possible to specify a whole pipeline of sorts
within one parameter.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

| **-o, --out**\ =\ *FILE*
| Write the sorted dynamic succinct variation graph to this file. A file
  ending with *.og* is recommended.

Files IO Options
---------------

| **-X, --path-index**\ =\ *FILE*
| Load the succinct variation graph index from this *FILE*. The file name usually ends with *.xp*.

| **-s, --sort-order**\ =\ *FILE*
| *FILE* containing the sort order. Each line contains one node
  identifier.

| **-C, --temp-dir**\ =\ *PATH*
| Directory for temporary files.

Topological Sort Options
-----------------

| **-b, --breadth-first**
| Use a (chunked) breadth first topological sort.

| **-B, --breadth-first-chunk**\ =\ *N*
| Chunk size for breadth first topological sort. Specify how many
  nucleotides to grap at once in each BFS phase.

| **-c, --cycle-breaking**
| Use a cycle breaking sort.

| **-z, --depth-first**
| Use a (chunked) depth first topological sort.

| **-Z, --depth-first-chunk**\ =\ *N*
| Chunk size for the depth first topological sort. Specify how many
  nucleotides to grap at once in each DFS phase.

| **-w, --two-way**
| Use a two-way topological algorithm for sorting. It is a maximum of
  head-first and tail-first topological sort.

| **-n, --no-seeds**
| Don’t use heads or tails to seed topological sort.

Random Sort Options
-----------

| **-r, --random**
| Randomly sort the graph.

DAGify Sort Options
-----------

| **-d, --dagify-sort**
| Sort on the basis of a DAGified graph.

Path Guided 1D Linear SGD Sort
------------------------------

| **-Y, --path-sgd**
| Apply the path-guided 1D linear SGD algorithm to organize the graph.

| **-f, --path-sgd-use-paths**\ =FILE
| Specify a line separated list of paths to sample from for the on the
  fly term generation process in the path guided linear 1D SGD (default: sample from all paths).

| **-G, --path-sgd-min-term-updates-paths**\ =\ *N*
| The minimum number of terms to be updated before a new path guided
  linear 1D SGD iteration with adjusted learning rate eta starts,
  expressed as a multiple of total path steps (default: *1.0*). Can be overwritten by *-U, -path-sgd-min-term-updates-nodes=N*.

| **-U, --path-sgd-min-term-updates-nodes**\ =\ *N*
| The minimum number of terms to be updated before a new path guided
  linear 1D SGD iteration with adjusted learning rate eta starts,
  expressed as a multiple of the number of nodes (default: NONE. *-G,path-sgd-min-term-updates-paths=N* is used).

| **-j, --path-sgd-delta**\ =\ *N*
| The threshold of maximum displacement approximately in bp at which to
  stop path guided linear 1D SGD (default: *0.0*).

| **-g, --path-sgd-eps**\ =\ *N*
| The final learning rate for path guided linear 1D SGD model (default: *0.01*).

| **-v, --path-sgd-eta-max**\ =\ *N*
| The first and maximum learning rate for path guided linear 1D SGD
  model (default: *squared steps of longest path in graph*).

| **-a, --path-sgd-zipf-theta**\ =\ *N*
| The theta value for the Zipfian distribution which is used as the
  sampling method for the second node of one term in the path guided
  linear 1D SGD model (default: *0.99*).

| **-x, --path-sgd-iter-max**\ =\ *N*
| The maximum number of iterations for path guided linear 1D SGD model (default: 30).

| **-F, --iteration-max-learning-rate**\ =\ *N*
| The iteration where the learning rate is max for path guided linear 1D SGD model (default: *0*).

| **-k, --path-sgd-zipf-space**\ =\ *N*
| The maximum space size of the Zipfian distribution which is used as
  the sampling method for the second node of one term in the path guided
  linear 1D SGD model (default: *longest path length*).

| **-I, --path-sgd-zipf-space-max**\ =\ *N*
| The maximum space size of the Zipfian distribution beyond which
  quantization occurs (default: *100*).

| **-l, --path-sgd-zipf-space-quantization-step**\ =\ *N*
| Quantization step size when the maximum space size of the Zipfian
  distribution is exceeded (default: *100*).

| **-y, --path-sgd-zipf-max-num-distributions**\ =\ *N*
| Approximate maximum number of Zipfian distributions to calculate (default: *100*).

| **-q, --path-sgd-seed**\ =\ *N*
| Set the seed for the deterministic 1-threaded path guided linear 1D SGD model (default: *pangenomic!*).

| **-u, --path-sgd-snapshot**\ =\ *STRING*
| Set the prefix to which each snapshot graph of a path guided 1D SGD
  iteration should be written to. This is turned off per default. This
  argument only works when *-Y, --path-sgd* was specified. Not applicable
  in a pipeline of sorts.

| **-H, --target-paths**\ =\ *FILE*
| Read the paths that should be considered as target paths (references) from this *FILE*. PG-SGD will keep the nodes of the given paths fixed. A path's rank determines it's weight for decision making and is given by its position in the given *FILE*.


Pipeline Sorting Options
----------------

| **-p, --pipeline**\ =\ *STRING*
| Apply a series of sorts, based on single character command line
  arguments given to this command (default: NONE). *s*: Topolocigal sort, heads only. *n*: Topological sort, no heads, no tails. *d*: DAGify sort. *c*: Cycle breaking sort. *b*: Breadth first topological sort. *z*: Depth first topological sort. *w*: Two-way topological sort. *r*: Random sort. *Y*: PG-SGD 1D sort. *f*: Reverse order. *g*: Groom the graph. An example could be *Ygs*.

Path Sorting Options
--------------------

| **-L, --paths-min**
| Sort paths by their lowest contained node identifier.

| **-M, --paths-max**
| Sort paths by their highest contained node identifier.

| **-A, --paths-avg**
| Sort paths by their average contained node identifier.

| **-R, --paths-avg-rev**
| Sort paths in reverse by their average contained node identifier.

| **-D, --path-delim**\ =\ *path-delim*
| Sort paths in bins by their prefix up to this delimiter.

Optimize Options
---------------------

| **-O, --optimize**
| Use the MutableHandleGraph::optimize method to compact the node
  identifier space.

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for the parallel operations.

Processing Information
----------------------

| **-P, --progress**
| Write the current progress to stderr.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi sort**.

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
