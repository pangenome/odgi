.. _odgi priv:

#########
odgi priv
#########

Differentially private sampling of graph subpaths, which applies the exponential mechanism to randomly sample shared sub-haplotypes with a given ε, target coverage, and minimum length.

SYNOPSIS
========

**odgi priv** [**-i, --idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

**odgi priv** generates a graph containing subpaths from the original graph which have been sampled under differential privacy guarantees.
We apply the exponential mechanism to guide our sampling of path intervals.
Each sampling iteration begins at a randomly-selected node and orientation.
All the steps which overlap this node are grouped by the next node which they reach, and the next step is chosen using the exponential mechanism.
To parameterize the exponential mechanism, we require a *utility*, which is a socially-defined value that we take to be the frequency of each subpath.
We then choose the next step with random sampling weighted by the exponentiated utility of each subpath extension scaled by a factor, ε (epsilon), which quantifies the degree of variation we expect when applying randomized algorithms to two datasets which differ by a single individual.
The process continues until a user-specified target length is reached, which defaults to 10kbp.

A few key parameters define properties of the output graph.
Low ε implies stronger privacy guarantees, while higher ε may include more rare haplotypes, and can be adjusted with **-e, --epsilon**.
Sampling continues until we reach a target depth of coverage over the graph, which defaults to 1x but can be changed by **-d, --target-depth**.
A minimum haplotype frequency is essential: the inclusion of singleton haplotypes violates differential privacy, so we do not recommend setting **-c, --min-hap-freq** below 2 except for testing the path cover properties of the system.
The target haplotype length **-b, --bp-target** is the minimum length haplotype to emit.
As very long haplotypes are rarely shared, except in extremely large cohorts, setting this too long (e.g. 100kb in 100 humans) will tend to cause the algorithm to stall, as it repeats sampling steps until it finds long haplotypes that occur at least **-c** times.
For optimal pangenome coverage, we suggest setting **-b** lower than the default, such as to 1kbp, with the caveat that this will reduce the length of novel variation that may be sampled from the graph.

Note that the output of **odgi priv** is a graph that does *not* meet the differential privacy guarantees.
We find that it is important to maintain the original graph for the purposes of debugging.
However, the full node space of the graph is preserved, which can lead to information leakage.
Further processing steps to prune 0-depth nodes and unchop the graph are suggested **odgi prune -i x.og -o - -c 1  | odgi sort -i - -o - -O | odgi unchop -i - -o - | odgi sort -i - -p Ygs -o y.og**.
The strongest privacy guarantees can be achieved by publishing only FASTA sequences for the paths of the emitted graph via **odgi paths -f**.


OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

| **-o, --out**\ =\ *FILE*
| Write the graph with sub-paths sampled under differential privacy to this FILE (.og recommended).

Differential Privacy Mechanism
------------------------------

| **-e, --epsilon**\ =\ *ε*
| Epsilon (ε) for exponential mechanism. [default: 0.01]

| **-d, --target-depth**\ =\ *DEPTH*
| Sample until we have approximately this path depth over the graph. [default: 1]

| **-c, --min-hap-freq**\ =\ *N*
| Minimum frequency (count) of haplotype observation to emit. Singularities occur at -c 1, so we warn against its use. [default: 2]

| **-b, --bp-target**\ =\ *bp*
| Target sampled haplotype length. All long haplotypes tend to be rare, so setting this to lengths greater than the typical recombination block size will result in long runtimes and poor sampling of the graph. [default: 10000]

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for parallel operations.

Processing Information
----------------------

| **-P, --progress**
| Print information about the operations and the progress to stderr.

| **-W, --write-haps**
| Write each sampled haplotype to stdout.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi paths**.

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
