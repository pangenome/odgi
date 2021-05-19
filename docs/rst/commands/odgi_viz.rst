.. _odgi viz:

#########
odgi viz
#########

variation graph visualizations

SYNOPSIS
========

**odgi viz** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
[*OPTION*]…

DESCRIPTION
===========

The odgi viz(1) command can produce a linear, static visualization of an
odgi variation graph. It can aggregate the pangenome into bins and
directly renders a raster image. The binning level can be specified in
input or it is calculated from the target width of the PNG to emit. Can
be used to produce visualizations for gigabase scale pangenomes. For
more information about the binning process, please refer to :ref:`odgi bin`.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph to convert from. The file
  name usually ends with *.og*.

| **-o, –out**\ =\ *FILE*
| Write the visualization in PNG format to this file.

Visualization Options
---------------------

| **-x, –width**\ =\ *N*
| Set the width in pixels of the output image.

| **-y, –height**\ =\ *N*
| Set the height in pixels of the output image.

| **-P, –path-height**\ =\ *N*
| The height in pixels for a path.

| **-X, –path-x-padding**\ =\ *N*
| The padding in pixels on the x-axis for a path.

| **-R, –pack-paths**
| Pack all paths rather than displaying a single path per row.

| **-L, –link-path-pieces**\ =\ *FLOAT*
| Show thin links of this relative width to connect path pieces.

| **-A, –alignment-prefix**\ =\ *STRING*
| Apply alignment related visual motifs to paths which have this name
  prefix. It affects the [**-S, –show-strand**] and [**-d,
  –change-darkness**] options.

| **-S, –show-strand**
| Use red and blue coloring to display forward and reverse alignments.
  This parameter can be set in combination with [**-A,
  –alignment-prefix**\ =\ *STRING*].

| **-z, –color-by-mean-inversion-rate**
| Change the color respect to the node strandness (black for forward,
  red for reverse); in binned mode (**-b, –binned-mode**), change the
  color respect to the mean inversion rate of the path for each bin,
  from black (no inversions) to red (bin mean inversion rate equals to
  1).

| **-s, –color-by-prefix**
| Colors paths by their names looking at the prefix before the given
  character C.

Intervals selection
-------------------

| **-r, –path-range**
| Nucleotide range to visualize: ``STRING=[PATH:]start-end``. ``\*-end``
  for ``[0,end]``; ``start-*`` for ``[start,pangenome_length]``. If no
  PATH is specified, the nucleotide positions refer to the pangenome’s
  sequence (i.e., the sequence obtained arranging all the graph’s node
  from left to right).

Paths selection
---------------

| **-p, –paths-to-display**
| List of paths to display in the specified order; the file must contain
  one path name per line and a subset of all paths can be specified.

Path names visualization Options
--------------------------------

| **-H, –hide-path-names**
| Hide the path names on the left of the generated image.

| **-C, –color-path-names-background**
| Color path names background with the same color as paths

| **-c, –max-num-of-characters**
| Maximum number of characters to display for each path name (max 128
  characters). The default value is *the length of the longest path
  name* (up to 32 characters).

Binned Mode Options
-------------------

| **-b, –binned-mode**
| The variation graph is binned before its visualization. Each pixel in
  the output image will correspond to a bin. For more information about
  the binning process, please refer to `odgi
  bin <#odgi_bin.adoc#_odgi_bin1>`__.

| **-w, –bin-width**\ =\ *N*
| The bin width specifies the size of each bin in the binned mode. If it
  is not specified, the bin width is calculated from the width in pixels
  of the output image.

| **-g, –no-gap-links**
| We divide links into 2 classes:

1. the links which help to follow complex variations. They need to be
   drawn, else one could not follow the sequence of a path.

2. the links helping to follow simple variations. These links are called
   **gap-links**. Such links solely connecting a path from left to right
   may not be relevant to understand a path’s traversal through the
   bins. Therefore, when this option is set, the gap-links are not drawn
   in binned mode.

| **-m, –color-by-mean-coverage**
| Change the color respect to the mean coverage of the path for each
  bin, from black (no coverage) to blue (max bin mean coverage in the
  entire graph).

Gradient Mode (also known as Position Mode) Options
---------------------------------------------------

| **-d, –change-darkness**
| Change the color darkness based on nucleotide position in the path.
  When it is used in binned mode, the mean inversion rate of the bin
  node is considered to set the color gradient starting position: when
  this rate is greater than 0.5, the bin is considered inverted, and the
  color gradient starts from the right-end of the bin. This parameter
  can be set in combination with [**-A,
  –alignment-prefix**\ =\ *STRING*].

| **-l, –longest-path**
| Use the longest path length to change the color darkness.

| **-u, –white-to-black**
| Change the color darkness from white (for the first nucleotide
  position) to black (for the last nucleotide position).

Program Information
-------------------

| **-h, –help**
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
