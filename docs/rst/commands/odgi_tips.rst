.. _odgi tips:

#########
odgi tips
#########

Identifying break point positions relative to given query (reference) path(s) of all the tips in the graph or of tips of given path(s). Prints BED records to stdout.

SYNOPSIS
========

**odgi tips** [**-i, --idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi tips command identifies break point positions relative to given target (reference) path(s) of all the tips in
the graph or of tips of given query path(s). Prints BED records to stdout. Each record consist of:

- **chrom**: The target path name.
- **start**: The 0-based start position of the query we hit in the node.
- **end**: The 1-based end position of the query we hit in the node.
- **path**: The name of the query path we walked.
- **path_pos**: The 0-based position of the query path we walked when we hit the node of the target path.
- **jaccard**: The jaccard index of the query and target path around the region of the step where the query hit the target.
- **walk_from_front**: If `1` we walked from the head of the target path. Else we walked from the tail and it is `0`.
- **additional_jaccards**: The additional jaccards of candidate reference step(s). Comma-separated.

OPTIONS
=======

MANDATORY OPTIONS
-----------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!


Tips Options
-------------

| **-q, --query-path**\ =\ *NAME*
| Use this query path.

| **-r, --target-path**\ =\ *NAME*
| Use this target (reference) path.

| **-Q, --query-paths**\ =\ *FILE*
| Use query paths listed in *FILE*. Each line must contain one path.

| **-R, --target-paths**\ =\ *FILE*
| Use target (reference) paths list (one per line) in *FILE*.

| **-v, --not-visited-tsv**\ =\ *FILE*
| Write query path(s) that do not visit the target path(s) to this *FILE*.

| **-n, --n-best**\ =\ *N*
| Report up to **N**th best target (reference) matches for each query path (default: 1).

| **-w, --jaccard-context**\ =\ *N*
| Maximum walking distance in nucleotides for one orientation when finding the best target (reference) range for each query path (default: 10000). Note: If we walked 9999 base pairs and **w, --jaccard-context** is **10000**, we will also include the next node, even if we overflow the actual limit.

| **-j, --jaccards**
| If for a target (reference) path several matches are possible, also report the additional jaccard indices (default: false). In the resulting BED, an '.' is added, if set to 'false'.

Step Index Options
------------------

| **-a, --step-index**\ =\ *FILE*
| Load the step index from this *FILE*. The file name usually ends with *.stpidx*. (default: build the step index from scratch with a sampling rate of 8).

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
