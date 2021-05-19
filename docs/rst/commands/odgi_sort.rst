.. _odgi sort:

#########
odgi sort
#########

sort a variation graph

SYNOPSIS
========

**odgi sort** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi sort(1) command sorts a succinct variation graph. Odgi sort
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

-  A 1D linear SGD sort: Odgi implements a 1D linear, variation graph
   adjusted, multi-threaded version of the `Graph Drawing by Stochastic
   Gradient Descent <https://arxiv.org/abs/1710.04626>`__ algorithm. The
   force-directed graph drawing algorithm minimizes the graph’s energy
   function or stress level. It applies stochastic gradient descent
   (SGD) to move a single pair of nodes at a time.

-  A path guided, 1D linear SGD sort: Odgi implements a 1D linear,
   variation graph adjusted, multi-threaded version of the `Graph
   Drawing by Stochastic Gradient
   Descent <https://arxiv.org/abs/1710.04626>`__ algorithm. The
   force-directed graph drawing algorithm minimizes the graph’s energy
   function or stress level. It applies stochastic gradient descent
   (SGD) to move a single pair of nodes at a time. The path index is
   used to pick the terms to move stochastically. If ran with 1 thread
   only, the resulting order of the graph is deterministic. The seed is
   adjustable.

Sorting the paths in a graph my refine the sorting process. For the
users’ convenience, it is possible to specify a whole pipeline of sorts
within one parameter.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to sort. The file name
  usually ends with *.og*.

| **-o, –out**\ =\ *FILE*
| Write the sorted dynamic succinct variation graph to this file. A file
  ending with *.og* is recommended.

| **-s, –sort-order**\ =\ *FILE*
| File containing the sort order. Each line contains one node
  identifier.

Topological Sorts
-----------------

| **-b, –breadth-first**
| Use a (chunked) breadth first topological sort.

| **-B, –breadth-first-chunk**\ =\ *N*
| Chunk size for breadth first topological sort. Specify how many
  nucleotides to grap at once in each BFS phase.

| **-z, –depth-first**
| Use a (chunked) depth first topological sort.

| **-Z, –depth-first-chunk**\ =\ *N*
| Chunk size for the depth first topological sort. Specify how many
  nucleotides to grap at once in each DFS phace.

| **-w, –two-way**
| Use a two-way topological algorithm for sorting. It is a maximum of
  head-first and tail-first topological sort.

| **-n, –no-seeds**
| Don’t use heads or tails to seed topological sort.

| **-c, –cycle-breaking**
| Use a cycle breaking sort.

Random Sort
-----------

| **-r, –random**
| Randomly sort the graph.

Path Guided 1D Linear SGD Sort
------------------------------

| **-Y, –path-sgd**
| Apply path guided 1D linear SGD algorithm to organize the graph.

| **-X, –path-index**\ =\ *FILE*
| Load the path index from this *FILE*.

| **-f, –path-sgd-use-paths**\ =FILE
| Specify a line separated list of paths to sample from for the on the
  fly term generation process in the path guided linear 1D SGD. The
  default value are *all paths*.

| **-G, –path-sgd-min-term-updates-paths**\ =\ *N*
| The minimum number of terms to be updated before a new path guided
  linear 1D SGD iteration with adjusted learning rate eta starts,
  expressed as a multiple of total path steps. The default value is
  *1.0*. Can be overwritten by *-U, -path-sgd-min-term-updates-nodes=N*.

| **-U, –path-sgd-min-term-updates-nodes**\ =\ *N*
| The minimum number of terms to be updated before a new path guided
  linear 1D SGD iteration with adjusted learning rate eta starts,
  expressed as a multiple of the number of nodes. Per default, the
  argument is not set. The default of *-G,
  path-sgd-min-term-updates-paths=N* is used).

| **-j, –path-sgd-delta**\ =\ *N*
| The threshold of maximum displacement approximately in bp at which to
  stop path guided linear 1D SGD. Default values is *0.0*.

| **-g, –path-sgd-eps**\ =\ *N*
| The final learning rate for path guided linear 1D SGD model. The
  default value is *0.01*.

| **-v, –path-sgd-eta-max**\ =\ *N*
| The first and maximum learning rate for path guided linear 1D SGD
  model. The default value is *squared steps of longest path in graph*.

| **-a, –path-sgd-zipf-theta**\ =\ *N*
| The theta value for the Zipfian distribution which is used as the
  sampling method for the second node of one term in the path guided
  linear 1D SGD model. The default value is *0.99*.

| **-x, –path-sgd-iter-max**\ =\ *N*
| The maximum number of iterations for path guided linear 1D SGD model.
  The default value is *30*.

| **-F, –iteration-max-learning-rate**\ =\ *N*
| The iteration where the learning rate is max for path guided linear 1D
  SGD model. The default value is *0*.

| **-k, –path-sgd-zipf-space**\ =\ *N*
| The maximum space size of the Zipfian distribution which is used as
  the sampling method for the second node of one term in the path guided
  linear 1D SGD model. The default value is the *longest path length*.

| **-I, –path-sgd-zipf-space-max**\ =\ *N*
| The maximum space size of the Zipfian distribution beyond which
  quantization occurs. Default value is *100*.

| **-l, –path-sgd-zipf-space-quantization-step**\ =\ *N*
| Quantization step size when the maximum space size of the Zipfian
  distribution is exceeded. Default value is *100*.

| **-y, –path-sgd-zipf-max-num-distributions**\ =\ *N*
| Approximate maximum number of Zipfian distributions to calculate. The
  default value is *100*.

| **-q, –path-sgd-seed**\ =\ *N*
| Set the seed for the deterministic 1-threaded path guided linear 1D
  SGD model. The default value is *pangenomic!*.

| **-u, –path-sgd-snapshot**\ =\ *STRING*
| Set the prefix to which each snapshot graph of a path guided 1D SGD
  iteration should be written to. This is turned off per default. This
  argument only works when *-Y, –path-sgd* was specified. Not applicable
  in a pipeline of sorts.

Path Sorting Options
--------------------

| **-L, –paths-min**
| Sort paths by their lowest contained node identifier.

| **-M, –paths-max**
| Sort paths by their highest contained node identifier.

| **-A, –paths-avg**
| Sort paths by their average contained node identifier.

| **-R, –paths-avg-rev**
| Sort paths in reverse by their average contained node identifier.

| **-D, –path-delim**\ =\ *path-delim*
| Sort paths in bins by their prefix up to this delimiter.

Pipeline Sorting
----------------

| **-p, –pipeline**\ =\ *STRING*
| Apply a series of sorts, based on single character command line
  arguments given to this command. The default sort is *s*. The reverse
  sort would be specified via *f*.

Additional Parameters
---------------------

| **-d, –dagify-sort**
| Sort on the basis of a DAGified graph.

| **-O, –Optimize**
| Use the MutableHandleGraph::optimize method to compact the node
  identifier space.

Threading
---------

| **-t, –threads**\ =\ *N*
| Number of threads to use for the parallel operations.

Processing Information
----------------------

| **-P, –progress**
| Print sort progress to stdout.

Program Information
-------------------

| **-h, –help**
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
