.. _odgi viz:

#########
odgi viz
#########

Visualize a variation graph in 1D.

SYNOPSIS
========

**odgi viz** [**-i, --idx**\ =\ *FILE*] [**-o, --out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi viz command can produce a linear, static visualization of an
odgi variation graph. It can aggregate the pangenome into bins and
directly renders a raster image. The binning level can be specified in
input or it is calculated from the target width of the PNG to emit. Can
be used to produce visualizations for gigabase scale pangenomes. For
more information about the binning process, please refer to :ref:`odgi bin`.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

| **-o, --out**\ =\ *FILE*
| Write the visualization in PNG format to this *FILE*.

Visualization Options
---------------------

| **-x, --width**\ =\ *N*
| Set the width in pixels of the output image (default: 1500).

| **-y, --height**\ =\ *N*
| Set the height in pixels of the output image (default: 500).

| **-a, --path-height**\ =\ *N*
| The height in pixels for a path.

| **-X, --path-x-padding**\ =\ *N*
| The padding in pixels on the x-axis for a path.

| **-n, --no-path-borders**
| Don't show path borders.

| **-b, --black-path-borders**
| Draw path borders in black (default: draw path borders in white).

| **-R, --pack-paths**
| Pack all paths rather than displaying a single path per row.

| **-L, --link-path-pieces**\ =\ *FLOAT*
| Show thin links of this relative width to connect path pieces.

| **-A, --alignment-prefix**\ =\ *STRING*
| Apply alignment related visual motifs to paths which have this name
  prefix. It affects the -S, --show-strand and -d, –change-darkness options.

| **-S, --show-strand**
| Use red and blue coloring to display forward and reverse alignments.
  This parameter can be set in combination with -A, –alignment-prefix=STRING.

| **-z, --color-by-mean-inversion-rate**
| Change the color respect to the node strandness (black for forward,
  red for reverse); in binned mode (**-b, --binned-mode**), change the
  color respect to the mean inversion rate of the path for each bin,
  from black (no inversions) to red (bin mean inversion rate equals to
  1).

| **-N, --color-by-uncalled-bases**
| Change the color with respect to the uncalled bases of the path for each
  bin, from black (no uncalled bases) to green (all uncalled bases).

| **-s, --color-by-prefix**\ =\ *CHAR*
| Color paths by their names looking at the prefix before the given
  character *CHAR*.

| **-M, --prefix-merges**\ =\ *FILE*
| Merge paths beginning with prefixes listed (one per line) in *FILE*.

| **-I, --ignore-prefix**\ =\ *PREFIX*
| Ignore paths starting with the given *PREFIX*.

Intervals Selection Options
-------------------

| **-r, --path-range**
| Nucleotide range to visualize: ``STRING=[PATH:]start-end``. ``*-end``
  for ``[0,end]``; ``start-*`` for ``[start,pangenome_length]``. If no
  PATH is specified, the nucleotide positions refer to the pangenome’s
  sequence (i.e., the sequence obtained arranging all the graph’s node
  from left to right).

Path Selection Options
---------------

| **-p, --paths-to-display**
| List of paths to display in the specified order; the file must contain
  one path name per line and a subset of all paths can be specified.

Path Names Viz Options
--------------------------------

| **-H, --hide-path-names**
| Hide the path names on the left of the generated image.

| **-C, --color-path-names-background**
| Color path names background with the same color as paths.

| **-c, --max-num-of-characters**\ =\ *N*
| Maximum number of characters to display for each path name (max 128
  characters). The default value is *the length of the longest path
  name* (up to 32 characters).

Binned Mode Options
-------------------

| **-w, --bin-width**\ =\ *N*
| The bin width specifies the size of each bin in the binned mode. If it
  is not specified, the bin width is calculated from the width in pixels
  of the output image.

| **-m, --color-by-mean-depth**
| Change the color with respect to the mean coverage of the path for each
  bin, from black (no coverage) to blue (max bin mean coverage in the
  entire graph).

| **-B, --colorbrewer-palette**\ =\ *SCHEME:N*
| Use the colorbrewer palette specified by the given *SCHEME*, with the number of levels *N*. Specifiy 'show' to see available palettes.

| **-G, --no-grey-depth**
| Use the colorbrewer palette specified for < 0.5x and ~1x coverage bins (default: these bins are light and neutral grey).

Gradient Mode Options
---------------------

| **-d, --change-darkness**
| Change the color darkness based on nucleotide position in the path.
  When it is used in binned mode, the mean inversion rate of the bin
  node is considered to set the color gradient starting position: when
  this rate is greater than 0.5, the bin is considered inverted, and the
  color gradient starts from the right-end of the bin. This parameter
  can be set in combination with -A, –alignment-prefix=*STRING*].

| **-l, --longest-path**
| Use the longest path length to change the color darkness.

| **-u, --white-to-black**
| Change the color darkness from white (for the first nucleotide
  position) to black (for the last nucleotide position).

Compressed Mode Options
-----------------------

| **-O, --compressed-mode**
| Compress the view vertically, summarizing the path coverage across all
  paths displaying the information using only one path **COMPRESSED_MODE**.
  A heatmap color-coding from https://colorbrewer2.org/#type=diverging&scheme=RdBu&n=11
  is used. Alternatively, one can enter a colorbrewer palette via -B, --colorbrewer-palette\ =\ *SCHEME:N*.

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
| Print a help message for **odgi viz**.

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
	Refer to the *odgi* issue tracker at https://github.com/pangenome/odgi/issues.
